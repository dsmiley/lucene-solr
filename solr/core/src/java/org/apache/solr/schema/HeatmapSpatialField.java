package org.apache.solr.schema;

import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

import org.apache.lucene.document.Field;
import org.apache.lucene.index.IndexOptions;
import org.apache.lucene.index.IndexReaderContext;
import org.apache.lucene.queries.function.ValueSource;
import org.apache.lucene.search.Query;
import org.apache.lucene.spatial.prefix.HeatmapFacetCounter;
import org.apache.lucene.spatial.prefix.IntersectsPrefixTreeQuery;
import org.apache.lucene.spatial.prefix.PrefixTreeFacetCounter;
import org.apache.lucene.spatial.prefix.PrefixTreeStrategy;
import org.apache.lucene.spatial.prefix.tree.Cell;
import org.apache.lucene.spatial.prefix.tree.CellIterator;
import org.apache.lucene.spatial.prefix.tree.SpatialPrefixTree;
import org.apache.lucene.spatial.query.SpatialArgs;
import org.apache.lucene.spatial.query.SpatialOperation;
import org.apache.lucene.spatial.query.UnsupportedSpatialOperation;
import org.apache.lucene.util.Bits;
import org.apache.lucene.util.BytesRef;
import org.apache.lucene.util.FixedBitSet;
import org.apache.solr.common.SolrException;
import org.apache.solr.util.MapListener;
import org.locationtech.spatial4j.shape.Point;
import org.locationtech.spatial4j.shape.Rectangle;
import org.locationtech.spatial4j.shape.Shape;

/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * A spatial field that is only for heatmap faceting.  It is much faster and customizable for this purpose than
 * {@link org.apache.lucene.spatial.prefix.PrefixTreeStrategy#calcFacets(IndexReaderContext, Bits, Shape, int, int)}.
 */
public class HeatmapSpatialField extends AbstractSpatialPrefixTreeFieldType<HeatmapSpatialField.HeatmapOnlyStrategy> {

  @Override
  protected void init(IndexSchema schema, Map<String, String> args) {
    Map<String, String> argsWithDefaults = new HashMap<>();
    //nocommit reference constants
    argsWithDefaults.put("geo", "true");
    argsWithDefaults.put("prefixTree", "packedQuad");
    argsWithDefaults.put("square", "true");
    argsWithDefaults.put("distanceUnits", "kilometers"); // redundant given geo=true but being clear is good
    argsWithDefaults.put("maxDistErr", "0.005");// rounds to nearest grid: 4.77m (level 23)
    argsWithDefaults.putAll(args);
    MapListener<String, String> argsWrap = new MapListener<>(argsWithDefaults);
    super.init(schema, argsWithDefaults);
    args.keySet().removeAll(argsWrap.getSeenKeys());

    String minDistErrStr = args.remove("minDistErr");
    if (minDistErrStr == null) {
      // 64 grid cells wide at world resolution: 360 / 64 * 111.20 = 626km
      minDistErrStr = "630";// rounds to nearest grid: 626km  (level 6)
    }
    double minDistErrDEG = distanceUnits.multiplierFromThisUnitToDegrees() * Double.parseDouble(minDistErrStr);
    int firstLevel = grid.getLevelForDistance(minDistErrDEG);
    usedLevelBits.clear(0, firstLevel); // end exclusive
    if (usedLevelBits.scanIsEmpty()) {
      throw new SolrException(SolrException.ErrorCode.BAD_REQUEST,
          "no active grid levels (bad minDistErr/maxDistErr?)");
    }

    // TODO offer alternative means to configure sparsely
    log.debug("prefixTree={} firstLevel={}", grid, usedLevelBits.nextSetBit(0));
  }

  @Override
  protected HeatmapOnlyStrategy newPrefixTreeStrategy(String fieldName) {
    return new HeatmapOnlyStrategy(grid, usedLevelBits, fieldName);
  }


  /** @lucene.internal */
  public static class HeatmapOnlyStrategy extends PrefixTreeStrategy {

    public static final org.apache.lucene.document.FieldType FIELD_TYPE = new org.apache.lucene.document.FieldType();
    static {
      FIELD_TYPE.setIndexOptions(IndexOptions.DOCS);
      FIELD_TYPE.setTokenized(false);// because at present we use a Field per token instead of only one Field for all tokens
      FIELD_TYPE.setOmitNorms(true);
      FIELD_TYPE.freeze();
    }

    private final String[] levelToFieldName; // sparse, only filled for levels we need
    private final int usedLevelCount;

    public HeatmapOnlyStrategy(SpatialPrefixTree grid, FixedBitSet usedLevels, String fieldName) {
      super(grid, fieldName);

      levelToFieldName = new String[grid.getMaxLevels() + 1]; // first slot is always empty
      usedLevelCount = usedLevels.cardinality();
      if (grid.getMaxLevels() + 1 != usedLevels.length()) {
        throw new IllegalArgumentException();
      }
      for (int level = 0; level <= grid.getMaxLevels(); level++) {
        if (usedLevels.get(level)) {
          levelToFieldName[level] = fieldName + "_L" + String.format(Locale.ROOT, "%02d", level);
        }
      }
      assert levelToFieldName[levelToFieldName.length - 1] != null : "last field shouldn't be empty";
    }

    public int levelOnOrAfter(int level) {
      // find next non-empty (last is always filled)
      while (levelToFieldName[level] == null) {
        level++;
      }
      return level;
    }

    @Override
    public Field[] createIndexableFields(Shape shape) {
      if (!(shape instanceof Point)) {
        throw new SolrException(SolrException.ErrorCode.BAD_REQUEST,
            getClass().getSimpleName() + " only supports indexing points; got: " + shape);
      }
      Field[] fields = new Field[usedLevelCount];
      int fieldsIdx = 0;
      // For each grid level, output a different Field
      CellIterator cellIterator = grid.getTreeCellIterator(shape, grid.getMaxLevels());
      Cell cell;
      do {
        cell = cellIterator.next();
        String fieldName = levelToFieldName[cell.getLevel()];
        if (fieldName != null) {
          BytesRef bytesRef = cell.getTokenBytesNoLeaf(null);//null: need new BytesRef
          assert fields[fieldsIdx] == null;
          fields[fieldsIdx++] = new Field(fieldName, bytesRef, FIELD_TYPE);
        }
      } while (cellIterator.hasNext());
      return fields;
    }

    @Override
    public ValueSource makeDistanceValueSource(Point queryPoint, double multiplier) {
      throw new UnsupportedOperationException("Not supported on " + HeatmapSpatialField.class.getSimpleName());
    }

    @Override
    public Query makeQuery(SpatialArgs args) {
      if (args.getOperation() != SpatialOperation.Intersects) {
        throw new UnsupportedSpatialOperation(args.getOperation());
      }
      Shape shape = args.getShape();
      int gridLevelIdx = levelOnOrAfter(grid.getLevelForDistance(args.resolveDistErr(ctx, SpatialArgs.DEFAULT_DISTERRPCT)));
      return makeIntersectsQueryAtLevel(shape, gridLevelIdx);
    }

    public Query makeIntersectsQueryAtLevel(Shape shape, int level) {
      boolean hasPrefixTerms = false;
      return new IntersectsPrefixTreeQuery(
          shape, getFieldNameForLevel(level), grid, level, level - 1, hasPrefixTerms);
    }

    public String getFieldNameForLevel(int level) {
      String fieldName = levelToFieldName[level];
      if (fieldName == null) {
        throw new IllegalArgumentException("bad facetLevel: " + level);
      }
      return fieldName;
    }

    /** Similar to {@link HeatmapFacetCounter#calcFacets(PrefixTreeStrategy, IndexReaderContext, Bits, Shape, int, int)}. */
    @Override
    public HeatmapFacetCounter.Heatmap calcFacets(IndexReaderContext context, Bits topAcceptDocs,
                                                  Shape inputShape, int facetLevel, int maxCells) throws IOException {
      String fieldName = getFieldNameForLevel(facetLevel);
      HeatmapFacetCounter.Heatmap heatmap = HeatmapFacetCounter.initHeatmap(inputShape, grid, facetLevel, maxCells);

      // see PrefixTreeFacetCounter ... but we do it differently
      final boolean hasPrefixTerms = false;
      PrefixTreeFacetCounter.compute(fieldName, grid,
          hasPrefixTerms, context, topAcceptDocs, inputShape, facetLevel,
          new PrefixTreeFacetCounter.FacetVisitor() {
            final double cellWidth = heatmap.region.getWidth() / heatmap.columns;
            final double cellHeight = heatmap.region.getHeight() / heatmap.rows;
            final double heatMinX = heatmap.region.getMinX();
            final double heatMinY = heatmap.region.getMinY();
            @Override
            public void visit(Cell cell, int count) {
              final Rectangle rect = (Rectangle) cell.getShape();
              if (cell.getLevel() == facetLevel) {//heatmap level; count it directly
                //convert to col & row
                int column;
                if (rect.getMinX() >= heatMinX) {
                  column = (int) Math.round((rect.getMinX() - heatMinX) / cellWidth);
                } else { // due to dateline wrap
                  column = (int) Math.round((rect.getMinX() + 360 - heatMinX) / cellWidth);
                }
                int row = (int) Math.round((rect.getMinY() - heatMinY) / cellHeight);
                //note: unfortunately, it's possible for us to visit adjacent cells to the heatmap (if the SpatialPrefixTree
                // allows adjacent cells to overlap on the seam), so we need to skip them
                if (column < 0 || column >= heatmap.columns || row < 0 || row >= heatmap.rows) {
                  return;
                }
                // increment
                heatmap.counts[column * heatmap.rows + row] += count;

              } else {
                throw new IllegalStateException("Unexpected & not supported yet.");
              }
            }
          });
      return heatmap;
    }
  } // class HeatmapOnlyStrategy

}
