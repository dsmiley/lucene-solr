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
package org.apache.solr.common.params;

import java.util.Locale;

public interface CollectionParams {
  /**
   * What action
   **/
  String ACTION = "action";
  String NAME = "name";


  enum LockLevel {
    CLUSTER(0),
    COLLECTION(1),
    SHARD(2),
    REPLICA(3),
    NONE(10);

    public final int level;

    LockLevel(int i) {
      this.level = i;
    }

    public LockLevel getChild() {
      return getLevel(level + 1);
    }

    public static LockLevel getLevel(int i) {
      for (LockLevel v : values()) {
        if (v.level == i) return v;
      }
      return null;
    }

    public boolean isHigherOrEqual(LockLevel that) {
      return that.level <= level;
    }
  }

  enum CollectionAction {
    CREATE(true, LockLevel.COLLECTION),
    DELETE(true, LockLevel.COLLECTION),
    RELOAD(true, LockLevel.COLLECTION),
    SYNCSHARD(true, LockLevel.SHARD),
    CREATEALIAS(true, LockLevel.COLLECTION),
    DELETEALIAS(true, LockLevel.COLLECTION),
    SPLITSHARD(true, LockLevel.SHARD),
    DELETESHARD(true, LockLevel.SHARD),
    CREATESHARD(true, LockLevel.COLLECTION),
    DELETEREPLICA(true, LockLevel.SHARD),
    FORCELEADER(true, LockLevel.SHARD),
    MIGRATE(true, LockLevel.SHARD),
    ADDROLE(true, LockLevel.NONE),
    REMOVEROLE(true, LockLevel.NONE),
    CLUSTERPROP(true, LockLevel.NONE),
    REQUESTSTATUS(false, LockLevel.NONE),
    DELETESTATUS(false, LockLevel.NONE),
    ADDREPLICA(true, LockLevel.SHARD),
    OVERSEERSTATUS(false, LockLevel.NONE),
    LIST(false, LockLevel.NONE),
    CLUSTERSTATUS(false, LockLevel.NONE),
    ADDREPLICAPROP(true, LockLevel.REPLICA),
    DELETEREPLICAPROP(true, LockLevel.REPLICA),
    BALANCESHARDUNIQUE(true, LockLevel.SHARD),
    REBALANCELEADERS(true, LockLevel.COLLECTION),
    MODIFYCOLLECTION(true, LockLevel.COLLECTION),
    MIGRATESTATEFORMAT(true, LockLevel.CLUSTER),
    //only for testing. it just waits for specified time
    // these are not exposed via collection API commands
    // but the overseer is aware of these tasks
    MOCK_COLL_TASK(false, LockLevel.COLLECTION),
    MOCK_SHARD_TASK(false, LockLevel.SHARD),
    MOCK_REPLICA_TASK(false, LockLevel.REPLICA)
    ;
    public final boolean isWrite;
    public final LockLevel lockLevel;

    CollectionAction(boolean isWrite, LockLevel level) {
      this.isWrite = isWrite;
      this.lockLevel = level;
    }

    public static CollectionAction get(String p) {
      if (p != null) {
        try {
          return CollectionAction.valueOf(p.toUpperCase(Locale.ROOT));
        } catch (Exception ex) {
        }
      }
      return null;
    }

    public boolean isEqual(String s) {
      if (s == null) return false;
      return toString().equals(s.toUpperCase(Locale.ROOT));
    }

    public String toLower() {
      return toString().toLowerCase(Locale.ROOT);
    }

  }
}
