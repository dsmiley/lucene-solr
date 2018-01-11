
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
package org.apache.solr.cloud;

import java.lang.invoke.MethodHandles;
import java.time.Instant;
import java.time.ZoneId;
import java.time.ZonedDateTime;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.time.temporal.ChronoUnit;
import java.time.temporal.TemporalAccessor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TimeZone;
import java.util.stream.Collectors;

import org.apache.solr.cloud.OverseerCollectionMessageHandler.Cmd;
import org.apache.solr.common.SolrException;
import org.apache.solr.common.cloud.ClusterState;
import org.apache.solr.common.cloud.ZkNodeProps;
import org.apache.solr.common.cloud.ZkStateReader;
import org.apache.solr.common.util.NamedList;
import org.apache.solr.common.util.StrUtils;
import org.apache.solr.update.processor.TimeRoutedAliasUpdateProcessor;
import org.apache.solr.util.DateMathParser;
import org.apache.solr.util.TimeZoneUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static java.time.format.DateTimeFormatter.ISO_INSTANT;
import static org.apache.solr.cloud.OverseerCollectionMessageHandler.COLL_CONF;
import static org.apache.solr.common.SolrException.ErrorCode.BAD_REQUEST;
import static org.apache.solr.common.params.CommonParams.NAME;
import static org.apache.solr.common.params.CommonParams.TZ;
import static org.apache.solr.handler.admin.CollectionsHandler.ROUTED_ALIAS_COLLECTION_PROP_PFX;
import static org.apache.solr.update.processor.TimeRoutedAliasUpdateProcessor.DATE_TIME_FORMATTER;


public class CreateAliasCmd implements Cmd {

  private static final Logger log = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  public static final String ROUTING_TYPE = "router.name";
  public static final String ROUTING_FIELD = "router.field";
  public static final String ROUTING_INCREMENT = "router.interval";
  public static final String ROUTING_MAX_FUTURE = "router.max-future-ms";
  public static final String START = "start";
  // Collection constants should all reflect names in the v2 structured input for this command, not v1
  // names used for CREATE
  public static final String CREATE_COLLECTION_CONFIG = "create-collection.config";
  public static final String CREATE_COLLECTION_ROUTER_NAME = "create-collection.router.name";
  public static final String CREATE_COLLECTION_ROUTER_FIELD = "create-collection.router.field";
  public static final String CREATE_COLLECTION_NUM_SHARDS = "create-collection.numShards";
  public static final String CREATE_COLLECTION_SHARDS = "create-collection.shards";
  public static final String CREATE_COLLECTION_REPLICATION_FACTOR = "create-collection.replicationFactor";
  public static final String CREATE_COLLECTION_NRT_REPLICAS = "create-collection.nrtReplicas";
  public static final String CREATE_COLLECTION_TLOG_REPLICAS = "create-collection.tlogReplicas";
  public static final String CREATE_COLLECTION_PULL_REPLICAS = "create-collection.pullReplicas";
  public static final String CREATE_COLLECTION_NODE_SET = "create-collection.nodeSet";
  public static final String CREATE_COLLECTION_SHUFFLE_NODES = "create-collection.shuffleNodes";
  public static final String CREATE_COLLECTION_MAX_SHARDS_PER_NODE = "create-collection.maxShardsPerNode";
  public static final String CREATE_COLLECTION_AUTO_ADD_REPLICAS = "create-collection.autoAddReplicas";
  public static final String CREATE_COLLECTION_RULE = "create-collection.rule";
  public static final String CREATE_COLLECTION_SNITCH = "create-collection.snitch";
  public static final String CREATE_COLLECTION_POLICY = "create-collection.policy";
  public static final String CREATE_COLLECTION_PROPERTIES = "create-collection.properties";


  private final OverseerCollectionMessageHandler ocmh;


  /**
   * Parameters required for creating a routed alias
   */
  public static final List<String> REQUIRED_ROUTING_PARAMS = Collections.unmodifiableList(Arrays.asList(
          NAME,
          START,
          ROUTING_FIELD,
          ROUTING_TYPE,
          ROUTING_INCREMENT));

  /**
   * Optional parameters for creating a routed alias excluding parameters for collection creation.
   */
  public static final List<String> NONREQUIRED_ROUTING_PARAMS = Collections.unmodifiableList(Arrays.asList(
          ROUTING_MAX_FUTURE,
          TZ));
  /**
   * Parameters used by routed Aliases to create collections.
   */
  public static final List<String> COLLECTION_ROUTING_PARAMS = Collections.unmodifiableList(Arrays.asList(
      CREATE_COLLECTION_CONFIG,
      CREATE_COLLECTION_ROUTER_FIELD,
      CREATE_COLLECTION_ROUTER_NAME,
      CREATE_COLLECTION_NUM_SHARDS,
      CREATE_COLLECTION_SHARDS,
      CREATE_COLLECTION_REPLICATION_FACTOR,
      CREATE_COLLECTION_NRT_REPLICAS,
      CREATE_COLLECTION_TLOG_REPLICAS,
      CREATE_COLLECTION_PULL_REPLICAS,
      CREATE_COLLECTION_NODE_SET,
      CREATE_COLLECTION_SHUFFLE_NODES,
      CREATE_COLLECTION_MAX_SHARDS_PER_NODE,
      CREATE_COLLECTION_AUTO_ADD_REPLICAS,
      CREATE_COLLECTION_RULE,
      CREATE_COLLECTION_SNITCH,
      CREATE_COLLECTION_POLICY,
      CREATE_COLLECTION_PROPERTIES));

  public static final List<String> CREATE_ROUTED_ALIAS_PARAMS;
  static {
    List<String> params = new ArrayList<>();
    params.addAll(REQUIRED_ROUTING_PARAMS);
    params.addAll(NONREQUIRED_ROUTING_PARAMS);
    params.addAll(COLLECTION_ROUTING_PARAMS);
    CREATE_ROUTED_ALIAS_PARAMS = Collections.unmodifiableList(params);
  }

  public CreateAliasCmd(OverseerCollectionMessageHandler ocmh) {
    this.ocmh = ocmh;
  }

  @Override
  public void call(ClusterState state, ZkNodeProps message, NamedList results)
      throws Exception {
    final String aliasName = message.getStr(NAME);
    ZkStateReader zkStateReader = ocmh.zkStateReader;
    ZkStateReader.AliasesManager holder = zkStateReader.aliasesHolder;
    if (!anyRoutingParams(message)) {
      final List<String> canonicalCollectionList = parseCollectionsParameter(message.get("collections"));
      final String canonicalCollectionsString = StrUtils.join(canonicalCollectionList, ',');
      validateAllCollectionsExistAndNoDups(canonicalCollectionList, zkStateReader);
      holder.applyModificationAndExportToZk(aliases -> aliases.cloneWithCollectionAlias(aliasName, canonicalCollectionsString));
    } else {
      final String routedField = message.getStr(ROUTING_FIELD);
      final String routingType = message.getStr(ROUTING_TYPE);
      final String tz = message.getStr(TZ);
      final String start = message.getStr(START);
      final String increment = message.getStr(ROUTING_INCREMENT);
      final String maxFutureMs = message.getStr(ROUTING_MAX_FUTURE);

      try {
        if (0 > Long.valueOf(maxFutureMs)) {
          throw new NumberFormatException("Negative value not allowed here");
        }
      } catch (NumberFormatException e) {
        throw new SolrException(BAD_REQUEST, ROUTING_MAX_FUTURE + " must be a valid long integer representing a number " +
            "of milliseconds greater than or equal to zero");
      }

      // Validate we got everything we need
      if (routedField == null || routingType == null || start == null || increment == null) {
        SolrException solrException = new SolrException(BAD_REQUEST, "If any of " + CREATE_ROUTED_ALIAS_PARAMS +
            " are supplied, then all of " + REQUIRED_ROUTING_PARAMS + " must be present.");
        log.error("Could not create routed alias",solrException);
        throw solrException;
      }

      if (!"time".equals(routingType)) {
        SolrException solrException = new SolrException(BAD_REQUEST, "Only time based routing is supported at this time");
        log.error("Could not create routed alias",solrException);
        throw solrException;
      }
      // Check for invalid timezone
      if(tz != null && !TimeZoneUtils.KNOWN_TIMEZONE_IDS.contains(tz)) {
        SolrException solrException = new SolrException(BAD_REQUEST, "Invalid timezone:" + tz);
        log.error("Could not create routed alias",solrException);
        throw solrException;

      }
      TimeZone zone;
      if (tz != null) {
        zone = TimeZoneUtils.getTimeZone(tz);
      } else {
        zone = TimeZoneUtils.getTimeZone("UTC");
      }
      DateTimeFormatter fmt = DATE_TIME_FORMATTER.withZone(zone.toZoneId());

      // check that the increment is valid date math
      String checkIncrement = ISO_INSTANT.format(Instant.now()) + increment;
      DateMathParser.parseMath(new Date(), checkIncrement); // exception if invalid increment

      Instant startTime = validateStart(zone, fmt, start);

      // check field
      String config = String.valueOf(message.getProperties().get(ROUTED_ALIAS_COLLECTION_PROP_PFX + COLL_CONF));
      if (!zkStateReader.getConfigManager().configExists(config)) {
        SolrException solrException = new SolrException(BAD_REQUEST, "Could not find config '" + config + "'");
        log.error("Could not create routed alias",solrException);
        throw solrException;
      }

      // It's too much work to check the routed field against the schema, there seems to be no good way to get
      // a copy of the schema aside from loading it directly from zookeeper based on the config name, but that
      // also requires I load solrconfig.xml to check what the value for managedSchemaResourceName is too, (or
      // discover that managed schema is not turned on and read schema.xml instead... and check for dynamic
      // field patterns too. As much as it would be nice to validate all inputs it's not worth the effort.

      String initialCollectionName = TimeRoutedAliasUpdateProcessor
          .formatCollectionNameFromInstant(aliasName, startTime, fmt );

      NamedList createResults = new NamedList();
      ZkNodeProps collectionProps = selectByPrefix(ROUTED_ALIAS_COLLECTION_PROP_PFX, message)
          .plus(NAME, initialCollectionName);
      Map<String, String> metadata = buildAliasMap(routedField, routingType, tz, increment, maxFutureMs, collectionProps);
      RoutedAliasCreateCollectionCmd.createCollectionAndWait(state,createResults,aliasName,metadata,initialCollectionName,ocmh);
      List<String> collectionList = Collections.singletonList(initialCollectionName);
      validateAllCollectionsExistAndNoDups(collectionList, zkStateReader);
      holder.applyModificationAndExportToZk(aliases -> aliases
          .cloneWithCollectionAlias(aliasName, initialCollectionName)
          .cloneWithCollectionAliasMetadata(aliasName, metadata));
    }

    // Sleep a bit to allow ZooKeeper state propagation.
    //
    // THIS IS A KLUDGE.
    //
    // Solr's view of the cluster is eventually consistent. *Eventually* all nodes and CloudSolrClients will be aware of
    // alias changes, but not immediately. If a newly created alias is queried, things should work right away since Solr
    // will attempt to see if it needs to get the latest aliases when it can't otherwise resolve the name.  However
    // modifications to an alias will take some time.
    //
    // We could levy this requirement on the client but they would probably always add an obligatory sleep, which is
    // just kicking the can down the road.  Perhaps ideally at this juncture here we could somehow wait until all
    // Solr nodes in the cluster have the latest aliases?
    Thread.sleep(100);
  }

  private Map<String, String> buildAliasMap(String routedField, String routingType, String tz, String increment, String maxFutureMs, ZkNodeProps collectionProps) {
    Map<String, Object> properties = collectionProps.getProperties();
    Map<String,String> cleanMap = properties.entrySet().stream()
        .filter(stringObjectEntry ->
            !"fromApi".equals(stringObjectEntry.getKey())
                && !"stateFormat".equals(stringObjectEntry.getKey())
                && !"name".equals(stringObjectEntry.getKey()))
        .collect(Collectors.toMap((e) -> "collection-create." + e.getKey(), e -> String.valueOf(e.getValue())));
    cleanMap.put(ROUTING_FIELD, routedField);
    cleanMap.put(ROUTING_TYPE, routingType);
    cleanMap.put(ROUTING_INCREMENT, increment);
    cleanMap.put(ROUTING_MAX_FUTURE, maxFutureMs);
    cleanMap.put(TZ, tz);
    return cleanMap;
  }

  private Instant validateStart(TimeZone zone, DateTimeFormatter fmt, String start) {
    // This is the normal/easy case, if we can get away with this great!
    TemporalAccessor startTime = attemptTimeStampParsing(start, zone.toZoneId());
    if (startTime == null) {
      // No luck, they gave us either date math, or garbage, so we have to do more work to figure out which and
      // to make sure it's valid date math and that it doesn't encode any millisecond precision.
      ZonedDateTime now = ZonedDateTime.now(zone.toZoneId()).truncatedTo(ChronoUnit.MILLIS);
      try {
        Date date = DateMathParser.parseMath(Date.from(now.toInstant()), start);
        String reformatted = fmt.format(date.toInstant().truncatedTo(ChronoUnit.MILLIS));
        Date reparse = Date.from(Instant.from(DATE_TIME_FORMATTER.parse(reformatted)));
        if (!reparse.equals(date)) {
          throw new SolrException(BAD_REQUEST,
              "Formatted time did not have the same milliseconds as original: " + date.getTime() + " vs. " +
                  reparse.getTime() + " This indicates that you used date math that includes milliseconds. " +
                  "(Hint: 'NOW' used without rounding always has this problem)" );
        }
        return date.toInstant();
      } catch (SolrException e) {
        throw new SolrException(BAD_REQUEST,
            "Start Time for the first collection must be a timestamp of the format yyyy-MM-dd_HH_mm_ss, " +
                "or a valid date math expression not containing specific milliseconds", e);
      }
    }
    return Instant.from(startTime);
  }

  private TemporalAccessor attemptTimeStampParsing(String start, ZoneId zone) {
    try {
      DATE_TIME_FORMATTER.withZone(zone);
      return DATE_TIME_FORMATTER.parse(start);
    } catch (DateTimeParseException e) {
      return null;
    }
  }

  private boolean anyRoutingParams(ZkNodeProps message) {

    return message.containsKey(ROUTING_FIELD) || message.containsKey(ROUTING_TYPE) || message.containsKey(START)
        || message.containsKey(ROUTING_INCREMENT) || message.containsKey(TZ);
  }

  private void validateAllCollectionsExistAndNoDups(List<String> collectionList, ZkStateReader zkStateReader) {
    final String collectionStr = StrUtils.join(collectionList, ',');

    if (new HashSet<>(collectionList).size() != collectionList.size()) {
      throw new SolrException(BAD_REQUEST,
          String.format(Locale.ROOT,  "Can't create collection alias for collections='%s', since it contains duplicates", collectionStr));
    }
    ClusterState clusterState = zkStateReader.getClusterState();
    Set<String> aliasNames = zkStateReader.getAliases().getCollectionAliasListMap().keySet();
    for (String collection : collectionList) {
      if (clusterState.getCollectionOrNull(collection) == null && !aliasNames.contains(collection)) {
        throw new SolrException(BAD_REQUEST,
            String.format(Locale.ROOT,  "Can't create collection alias for collections='%s', '%s' is not an existing collection or alias", collectionStr, collection));
      }
    }
  }

  /**
   * The v2 API directs that the 'collections' parameter be provided as a JSON array (e.g. ["a", "b"]).  We also
   * maintain support for the legacy format, a comma-separated list (e.g. a,b).
   */
  @SuppressWarnings("unchecked")
  private List<String> parseCollectionsParameter(Object colls) {
    if (colls == null) throw new SolrException(BAD_REQUEST, "missing collections param");
    if (colls instanceof List) return (List<String>) colls;
    return StrUtils.splitSmart(colls.toString(), ",", true).stream()
        .map(String::trim)
        .collect(Collectors.toList());
  }

  private ZkNodeProps selectByPrefix(String prefix, ZkNodeProps source) {
    final ZkNodeProps[] subSet = {new ZkNodeProps()};
    source.getProperties().entrySet().stream()
        .filter(entry -> entry.getKey().startsWith(prefix))
        .forEach(e -> subSet[0] = subSet[0].plus(e.getKey().substring(prefix.length()),e.getValue()));
    return subSet[0];
  }

}
