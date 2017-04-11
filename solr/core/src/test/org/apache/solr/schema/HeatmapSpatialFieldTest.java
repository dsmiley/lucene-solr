package org.apache.solr.schema;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.carrotsearch.randomizedtesting.annotations.Repeat;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.TotalHitCountCollector;
import org.apache.lucene.spatial.prefix.HeatmapFacetCounter;
import org.apache.lucene.util.Bits;
import org.apache.solr.SolrTestCaseJ4;
import org.apache.solr.search.SolrIndexSearcher;
import org.apache.solr.util.RefCounted;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.locationtech.spatial4j.context.SpatialContext;
import org.locationtech.spatial4j.distance.DistanceUtils;
import org.locationtech.spatial4j.shape.Point;
import org.locationtech.spatial4j.shape.Rectangle;
import org.locationtech.spatial4j.shape.Shape;
import org.locationtech.spatial4j.shape.SpatialRelation;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static com.carrotsearch.randomizedtesting.RandomizedTest.atMost;
import static com.carrotsearch.randomizedtesting.RandomizedTest.randomDouble;
import static com.carrotsearch.randomizedtesting.RandomizedTest.randomGaussian;
import static com.carrotsearch.randomizedtesting.RandomizedTest.randomInt;
import static com.carrotsearch.randomizedtesting.RandomizedTest.randomIntBetween;

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

public class HeatmapSpatialFieldTest extends SolrTestCaseJ4 {

  private static final Logger log = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private static File tmpSolrHome;
  private static File tmpConfDir;

  private static final String collection = "collection1";
  private static final String confDir = collection + "/conf";

  int cellsValidated;
  int cellValidatedNonZero;

  static final String ID_FLD = "str";
  static final String HM_FLD = "myHeatmapField";
  private static HeatmapSpatialField.HeatmapOnlyStrategy strategy;
  private static SpatialContext ctx;

  private RefCounted<SolrIndexSearcher> searcherRefCounted;

  /** TODO @see SpatialRPTFieldTypeTest */
  @BeforeClass
  private static void initManagedSchemaCore() throws Exception {
    // This testing approach means no new solrconfig or schema file or per-test temp solr-home!
    System.setProperty("managed.schema.mutable", "true");
    System.setProperty("managed.schema.resourceName", "schema-one-field-no-dynamic-field-unique-key.xml");
    System.setProperty("enable.update.log", "false");
    initCore("solrconfig-managed-schema.xml", "ignoredSchemaName?");

    IndexSchema oldSchema = h.getCore().getLatestSchema();

    HeatmapSpatialField fieldType = new HeatmapSpatialField();
    Map<String, String> ftInitMap = new HashMap<>();
    ftInitMap.put("prefixTree", "packedQuad");
    ftInitMap.put("square", "true");
    ftInitMap.put("minDistErr", "1000");
    ftInitMap.put("maxDistErr", "50");
    ftInitMap.put("distanceUnits", "kilometers");
    fieldType.init(oldSchema, ftInitMap);
    fieldType.setTypeName("heatmapType");
    SchemaField schemaField = new SchemaField(HM_FLD, fieldType, SchemaField.STORED | SchemaField.INDEXED, null);
    boolean persist = false; // don't write to test resource dir
    IndexSchema newSchema = oldSchema.addField(schemaField, persist);

    h.getCore().setLatestSchema(newSchema);

    strategy = fieldType.getStrategy(schemaField.getName());
    ctx = strategy.getSpatialContext();
  }

  private SolrIndexSearcher getIndexSearcher() {
    if (this.searcherRefCounted == null) {
      this.searcherRefCounted = h.getCore().getSearcher();
    }
    return this.searcherRefCounted.get();
  }

  private void closeIndexSearcher() {
    if (this.searcherRefCounted != null) {
      this.searcherRefCounted.decref();
      this.searcherRefCounted = null;
    }
  }

  @After
  public void after() {
    closeIndexSearcher();
    log.info("Validated " + cellsValidated + " cells, " + cellValidatedNonZero + " non-zero");
  }


  private void adoc(String id, Shape shape) {
    assertU(adoc(ID_FLD, id, HM_FLD, ctx.getFormats().getWktWriter().toString(shape)));
  }

  private void commit() {
    closeIndexSearcher();
    assertU(super.commit());
  }

  void deleteDoc(String id) {
    assertU(delI(id));
  }

//  @Test
//  public void testQueryCircle() throws IOException {
//    //overwrite setUp; non-geo bounds is more straight-forward; otherwise 88,88 would actually be practically north,
//    final SpatialContextFactory spatialContextFactory = new SpatialContextFactory();
//    spatialContextFactory.geo = false;
//    spatialContextFactory.worldBounds = new RectangleImpl(-90, 90, -90, 90, null);
//    ctx = spatialContextFactory.newSpatialContext();
//    final int LEVEL = 4;
//    grid = new QuadPrefixTree(ctx, LEVEL);
//    strategy = new RecursivePrefixTreeStrategy(grid, getTestClass().getSimpleName());
//    Circle circle = ctx.makeCircle(0, 0, 89);
//    adoc("0", ctx.makePoint(88, 88));//top-right, inside bbox of circle but not the circle
//    adoc("1", ctx.makePoint(0, 0));//clearly inside; dead center in fact
//    commit();
//    final HeatmapFacetCounter.Heatmap heatmap = HeatmapFacetCounter.calcFacets(
//        (PrefixTreeStrategy) strategy, indexSearcher.getTopReaderContext(), null,
//        circle, LEVEL, 1000);
//    //assert that only one point is found, not 2
//    boolean foundOne = false;
//    for (int count : heatmap.counts) {
//      switch (count) {
//        case 0: break;
//        case 1:
//          assertFalse(foundOne);//this is the first
//          foundOne = true;
//          break;
//        default:
//          fail("counts should be 0 or 1: " + count);
//      }
//    }
//    assertTrue(foundOne);
//  }

  @Test
  @Repeat(iterations = 20)
  public void testRandom() throws IOException {
    // Tests using random index shapes & query shapes. This has found all sorts of edge case bugs (e.g. dateline,
    // cell border, overflow(?)).

    final int numIndexedShapes = 1 + atMost(9);
    List<Shape> indexedShapes = new ArrayList<>(numIndexedShapes);
    for (int i = 0; i < numIndexedShapes; i++) {
      indexedShapes.add(randomIndexedShape());
    }

    //Main index loop:
    for (int i = 0; i < indexedShapes.size(); i++) {
      Shape shape = indexedShapes.get(i);
      adoc("" + i, shape);

      if (random().nextInt(10) == 0)
        commit();//intermediate commit, produces extra segments
    }
    //delete some documents randomly
    for (int id = 0; id < indexedShapes.size(); id++) {
      if (random().nextInt(10) == 0) {
        deleteDoc("" + id);
        indexedShapes.set(id, null);
      }
    }

    commit();

    // once without dateline wrap
    final Rectangle rect = randomRectangle();
    int firstLevel = strategy.levelOnOrAfter(1);
    queryHeatmapRecursive(usually() ? ctx.getWorldBounds() : rect, firstLevel);
    // and once with dateline wrap
    if (rect.getWidth() > 0) {
      double shift = random().nextDouble() % rect.getWidth();
      queryHeatmapRecursive(ctx.getShapeFactory().rect(
          DistanceUtils.normLonDEG(rect.getMinX() - shift),
          DistanceUtils.normLonDEG(rect.getMaxX() - shift),
          rect.getMinY(), rect.getMaxY()),
          firstLevel);
    }
  }

  /** Build heatmap, validate results, then descend recursively to another facet level. */
  private boolean queryHeatmapRecursive(Rectangle inputRange, int facetLevel) throws IOException {
    if (!inputRange.hasArea()) {
      // Don't test line inputs. It's not that we don't support it but it is more challenging to test if per-chance it
      // coincides with a grid line due due to edge overlap issue for some grid implementations (geo & quad).
      return false;
    }
    Bits filter = null; //FYI testing filtering of underlying PrefixTreeFacetCounter is done in another test
    //Calculate facets
    final int maxCells = 10_000;
    final HeatmapFacetCounter.Heatmap heatmap =
        strategy.calcFacets(getIndexSearcher().getTopReaderContext(), filter, inputRange, facetLevel, maxCells);

    validateHeatmapResult(inputRange, facetLevel, heatmap);

    boolean foundNonZeroCount = false;
    for (int count : heatmap.counts) {
      if (count > 0) {
        foundNonZeroCount = true;
        break;
      }
    }

    //Test again recursively to higher facetLevel (more detailed cells)
    if (foundNonZeroCount && cellsValidated <= 500 && facetLevel != strategy.getGrid().getMaxLevels() && inputRange.hasArea()) {
      for (int i = 0; i < 5; i++) {//try multiple times until we find non-zero counts
        if (queryHeatmapRecursive(randomRectangle(inputRange), strategy.levelOnOrAfter(facetLevel + 1))) {
          break;//we found data here so we needn't try again
        }
      }
    }
    return foundNonZeroCount;
  }

  private void validateHeatmapResult(Rectangle inputRange, int facetLevel, HeatmapFacetCounter.Heatmap heatmap)
      throws IOException {
    final Rectangle heatRect = heatmap.region;
    assertTrue(inputRange.relate(heatRect) == SpatialRelation.WITHIN || heatRect.equals(inputRange));
    final double cellWidth = heatRect.getWidth() / heatmap.columns;
    final double cellHeight = heatRect.getHeight() / heatmap.rows;
    for (int c = 0; c < heatmap.columns; c++) {
      for (int r = 0; r < heatmap.rows; r++) {
        final int facetCount = heatmap.getCount(c, r);
        double x = DistanceUtils.normLonDEG(heatRect.getMinX() + c * cellWidth + cellWidth / 2);
        double y = DistanceUtils.normLatDEG(heatRect.getMinY() + r * cellHeight + cellHeight / 2);
        Point pt = heatRect.getContext().getShapeFactory().pointXY(x, y);
        assertEquals("input:" + inputRange + " level:" + facetLevel, countMatchingDocsAtLevel(pt, facetLevel), facetCount);
      }
    }
  }

  private int countMatchingDocsAtLevel(Point pt, int facetLevel) throws IOException {
    Query filter = strategy.makeIntersectsQueryAtLevel(pt, facetLevel);
    TotalHitCountCollector collector = new TotalHitCountCollector();
    getIndexSearcher().search(filter, collector);
    cellsValidated++;
    if (collector.getTotalHits() > 0) {
      cellValidatedNonZero++;
    }
    return collector.getTotalHits();
  }

  private Shape randomIndexedShape() {
    return randomPoint();
  }

  //
  //  COPIED FROM SpatialTestCase:
  //

  protected Point randomPoint() {
    final Rectangle WB = ctx.getWorldBounds();
    return ctx.makePoint(
        randomIntBetween((int) WB.getMinX(), (int) WB.getMaxX()),
        randomIntBetween((int) WB.getMinY(), (int) WB.getMaxY()));
  }

  protected Rectangle randomRectangle() {
    return randomRectangle(ctx.getWorldBounds());
  }

  protected Rectangle randomRectangle(Rectangle bounds) {
    double[] xNewStartAndWidth = randomSubRange(bounds.getMinX(), bounds.getWidth());
    double xMin = xNewStartAndWidth[0];
    double xMax = xMin + xNewStartAndWidth[1];
    if (bounds.getCrossesDateLine()) {
      xMin = DistanceUtils.normLonDEG(xMin);
      xMax = DistanceUtils.normLonDEG(xMax);
    }

    double[] yNewStartAndHeight = randomSubRange(bounds.getMinY(), bounds.getHeight());
    double yMin = yNewStartAndHeight[0];
    double yMax = yMin + yNewStartAndHeight[1];

    return ctx.makeRectangle(xMin, xMax, yMin, yMax);
  }

  /** Returns new minStart and new length that is inside the range specified by the arguments. */
  protected double[] randomSubRange(double boundStart, double boundLen) {
    if (boundLen >= 3 && usually()) { // typical
      // prefer integers for ease of debugability ... and prefer 1/16th of bound
      int intBoundStart = (int) Math.ceil(boundStart);
      int intBoundEnd = (int) (boundStart + boundLen);
      int intBoundLen = intBoundEnd - intBoundStart;
      int newLen = (int) randomGaussianMeanMax(intBoundLen / 16.0, intBoundLen);
      int newStart = intBoundStart + randomInt(intBoundLen - newLen);
      return new double[]{newStart, newLen};
    } else { // (no int rounding)
      double newLen = randomGaussianMeanMax(boundLen / 16, boundLen);
      double newStart = boundStart + (boundLen - newLen == 0 ? 0 : (randomDouble() % (boundLen - newLen)));
      return new double[]{newStart, newLen};
    }
  }

  private double randomGaussianMinMeanMax(double min, double mean, double max) {
    assert mean > min;
    return randomGaussianMeanMax(mean - min, max - min) + min;
  }

  /**
   * Within one standard deviation (68% of the time) the result is "close" to
   * mean. By "close": when greater than mean, it's the lesser of 2*mean or half
   * way to max, when lesser than mean, it's the greater of max-2*mean or half
   * way to 0. The other 32% of the time it's in the rest of the range, touching
   * either 0 or max but never exceeding.
   */
  private double randomGaussianMeanMax(double mean, double max) {
    // DWS: I verified the results empirically
    assert mean <= max && mean >= 0;
    double g = randomGaussian();
    double mean2 = mean;
    double flip = 1;
    if (g < 0) {
      mean2 = max - mean;
      flip = -1;
      g *= -1;
    }
    // pivot is the distance from mean2 towards max where the boundary of
    // 1 standard deviation alters the calculation
    double pivotMax = max - mean2;
    double pivot = Math.min(mean2, pivotMax / 2);//from 0 to max-mean2
    assert pivot >= 0 && pivotMax >= pivot && g >= 0;
    double pivotResult;
    if (g <= 1)
      pivotResult = pivot * g;
    else
      pivotResult = Math.min(pivotMax, (g - 1) * (pivotMax - pivot) + pivot);

    double result = mean + flip * pivotResult;
    return (result < 0 || result > max) ? mean : result; // due this due to computational numerical precision
  }
}
