// === JMA2001走時表の実装 ===
class TravelTimeTable {
  constructor() {
    // tjma2001h.txtのデータを格納
    // フォーマット: P波走時(秒), S波走時(秒), 深さ(km), 震央距離(km)
    this.table = [];
    this.loaded = false;
  }

  // 走時表データをパース
  parseTableData(tableText) {
    const lines = tableText.trim().split('\n');
    this.table = [];

    for (const line of lines) {
      const parts = line.trim().split(/\s+/);
      if (parts.length >= 8 && parts[0] === 'P') {
        const pTime = parseFloat(parts[1]);
        const sTime = parseFloat(parts[3]);
        const depth = parseFloat(parts[4]);
        const distance = parseFloat(parts[5]);

        this.table.push({
          pTime,
          sTime,
          depth,
          distance
        });
      }
    }

    this.loaded = true;
    console.log(`走時表データ読み込み完了: ${this.table.length}エントリ`);
  }

  // 走時表データを設定(手動)
  setTableData(data) {
    this.table = data;
    this.loaded = true;
  }

  // 2点間の線形補間
  linearInterpolate(x, x0, x1, y0, y1) {
    if (x1 === x0) return y0;
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }

  // バイリニア補間(距離と深さの2次元補間)
  bilinearInterpolate(distance, depth, waveType = 'P') {
    if (!this.loaded || this.table.length === 0) {
      // フォールバック: 簡易計算
      const avgVelocity = waveType === 'P' ? 7.0 : 4.0;
      return distance / avgVelocity + depth / 8.0;
    }

    // 距離と深さの範囲を取得
    const distances = [...new Set(this.table.map(e => e.distance))].sort((a, b) => a - b);
    const depths = [...new Set(this.table.map(e => e.depth))].sort((a, b) => a - b);

    // 最も近い4点を探す
    let d0 = distances[0], d1 = distances[0];
    let z0 = depths[0], z1 = depths[0];

    // 距離の範囲を見つける
    for (let i = 0; i < distances.length - 1; i++) {
      if (distance >= distances[i] && distance <= distances[i + 1]) {
        d0 = distances[i];
        d1 = distances[i + 1];
        break;
      }
    }
    if (distance > distances[distances.length - 1]) {
      d0 = distances[distances.length - 2];
      d1 = distances[distances.length - 1];
    }

    // 深さの範囲を見つける
    for (let i = 0; i < depths.length - 1; i++) {
      if (depth >= depths[i] && depth <= depths[i + 1]) {
        z0 = depths[i];
        z1 = depths[i + 1];
        break;
      }
    }
    if (depth > depths[depths.length - 1]) {
      z0 = depths[depths.length - 2];
      z1 = depths[depths.length - 1];
    }

    // 4隅の点を取得
    const getPoint = (dist, dep) => {
      return this.table.find(e => e.distance === dist && e.depth === dep);
    };

    const p00 = getPoint(d0, z0);
    const p01 = getPoint(d0, z1);
    const p10 = getPoint(d1, z0);
    const p11 = getPoint(d1, z1);

    // いずれかの点が見つからない場合は最近傍法
    if (!p00 || !p01 || !p10 || !p11) {
      let nearest = this.table[0];
      let minDist = Infinity;

      for (const point of this.table) {
        const dist = Math.sqrt(
          Math.pow(point.distance - distance, 2) + 
          Math.pow(point.depth - depth, 2)
        );
        if (dist < minDist) {
          minDist = dist;
          nearest = point;
        }
      }

      return waveType === 'P' ? nearest.pTime : nearest.sTime;
    }

    // バイリニア補間
    const timeKey = waveType === 'P' ? 'pTime' : 'sTime';
    
    const t0 = this.linearInterpolate(distance, d0, d1, p00[timeKey], p10[timeKey]);
    const t1 = this.linearInterpolate(distance, d0, d1, p01[timeKey], p11[timeKey]);
    const result = this.linearInterpolate(depth, z0, z1, t0, t1);

    return result;
  }

  // P波の走時を取得
  getPWaveTravelTime(epicentralDistance, depth) {
    return this.bilinearInterpolate(epicentralDistance, depth, 'P');
  }

  // S波の走時を取得
  getSWaveTravelTime(epicentralDistance, depth) {
    return this.bilinearInterpolate(epicentralDistance, depth, 'S');
  }
}

// グローバル走時表インスタンス
const globalTravelTimeTable = new TravelTimeTable();

// tjma2001h.txtの基本データ(一部)を初期設定
// 実際のファイルを読み込む場合は parseTableData() を使用
globalTravelTimeTable.setTableData([
  { pTime: 0.082, sTime: 0.139, depth: 0, distance: 0 },
  { pTime: 0.420, sTime: 0.708, depth: 0, distance: 2 },
  { pTime: 0.825, sTime: 1.393, depth: 0, distance: 4 },
  { pTime: 1.231, sTime: 2.078, depth: 0, distance: 6 },
  { pTime: 1.633, sTime: 2.758, depth: 0, distance: 8 },
  { pTime: 2.030, sTime: 3.429, depth: 0, distance: 10 },
  { pTime: 2.422, sTime: 4.092, depth: 0, distance: 12 },
  { pTime: 2.809, sTime: 4.746, depth: 0, distance: 14 },
  // より多くのデータを追加する必要があります
  // 深さ別、距離別のデータが必要
]);

// === 震源推定クラス ===
class EpicenterEstimator {
  constructor(config = {}) {
    this.config = {
      initialDepth: 10, // 初期深さ(km)
      initialTimeOffset: -2, // 初期発生時刻オフセット(秒)
      gridStepLarge: 0.5, // 大まかな探索ステップ(度)
      gridStepSmall: 0.1, // 細かい探索ステップ(度)
      depthStepLarge: 50, // 大まかな深さステップ(km)
      depthStepSmall: 10, // 細かい深さステップ(km)
      nearStationWeight: 1.0, // 近傍観測点の重み
      nearStationRadius: 50, // 近傍観測点の半径(km)
      minDetections: 3, // 最小検知数
      undetectedPenalty: 1.0, // 未検知ペナルティ
      undetectedCheckTime: 3, // 未検知チェック時間(秒)
      ...config
    };
    
    this.currentEpicenter = null;
  }

  // 震央距離を計算(度から距離へ変換)
  calculateEpicentralDistance(lat1, lon1, lat2, lon2) {
    const R = 6371; // 地球の半径(km)
    const dLat = (lat2 - lat1) * Math.PI / 180;
    const dLon = (lon2 - lon1) * Math.PI / 180;
    const a = Math.sin(dLat/2) * Math.sin(dLat/2) +
              Math.cos(lat1 * Math.PI / 180) * Math.cos(lat2 * Math.PI / 180) *
              Math.sin(dLon/2) * Math.sin(dLon/2);
    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
    return R * c;
  }

  // 誤差レベルを計算
  calculateErrorLevel(epicenter, detections, currentTime) {
    const { lat, lon, depth, originTime } = epicenter;
    const originTimes = [];
    const weights = [];
    let errorLevel = 0;

    // 各観測点で発生時刻を計算
    detections.forEach(detection => {
      const distance = this.calculateEpicentralDistance(
        lat, lon,
        detection.location[0], detection.location[1]
      );
      
      const travelTime = globalTravelTimeTable.getPWaveTravelTime(distance, depth);
      const calculatedOriginTime = detection.detectionTime - travelTime;
      
      // 重みを計算(近い観測点ほど重い)
      let weight = 1.0;
      if (detections.length > 0) {
        const firstDistance = this.calculateEpicentralDistance(
          lat, lon,
          detections[0].location[0], detections[0].location[1]
        );
        weight = distance < this.config.nearStationRadius 
          ? 1.0 
          : firstDistance / distance;
      }
      
      originTimes.push(calculatedOriginTime);
      weights.push(weight);
    });

    // 発生時刻の平均を計算
    const weightSum = weights.reduce((a, b) => a + b, 0);
    const avgOriginTime = originTimes.reduce((sum, time, i) => 
      sum + time * weights[i], 0) / weightSum;

    // 2乗誤差を計算
    originTimes.forEach((time, i) => {
      const diff = time - avgOriginTime;
      errorLevel += (diff * diff) * weights[i];
    });

    // 着未着法による補正(揺れ検知初期のみ)
    const elapsedTime = (currentTime - detections[0].detectionTime) / 1000;
    if (elapsedTime < this.config.undetectedCheckTime || 
        detections.length < 30) {
      
      // 未検知観測点をチェック
      // (実装には全観測点リストが必要 - ここでは省略)
      // 計算上到達済みなのに未検知の観測点があればペナルティ
    }

    return errorLevel;
  }

  // 震源候補を生成
  generateCandidates(center, latStep, lonStep, depthStep = null) {
    const candidates = [];
    
    // 緯度・経度の候補
    candidates.push({
      lat: center.lat + latStep,
      lon: center.lon,
      depth: center.depth,
      originTime: center.originTime
    });
    candidates.push({
      lat: center.lat - latStep,
      lon: center.lon,
      depth: center.depth,
      originTime: center.originTime
    });
    candidates.push({
      lat: center.lat,
      lon: center.lon + lonStep,
      depth: center.depth,
      originTime: center.originTime
    });
    candidates.push({
      lat: center.lat,
      lon: center.lon - lonStep,
      depth: center.depth,
      originTime: center.originTime
    });

    // 深さの候補(指定された場合)
    if (depthStep !== null) {
      candidates.push({
        lat: center.lat,
        lon: center.lon,
        depth: center.depth + depthStep,
        originTime: center.originTime
      });
      candidates.push({
        lat: center.lat,
        lon: center.lon,
        depth: Math.max(0, center.depth - depthStep),
        originTime: center.originTime
      });
    }

    return candidates;
  }

  // 震源を反復的に改善
  improveEpicenter(initialEpicenter, detections, currentTime, 
                   latStep, lonStep, depthStep = null) {
    let currentEpicenter = { ...initialEpicenter };
    let improved = true;

    while (improved) {
      improved = false;
      const currentError = this.calculateErrorLevel(
        currentEpicenter, detections, currentTime
      );

      const candidates = this.generateCandidates(
        currentEpicenter, latStep, lonStep, depthStep
      );

      for (const candidate of candidates) {
        const candidateError = this.calculateErrorLevel(
          candidate, detections, currentTime
        );

        if (candidateError < currentError) {
          currentEpicenter = { ...candidate };
          improved = true;
          break;
        }
      }
    }

    return currentEpicenter;
  }

  // メイン震源推定処理
  estimate(detections, currentTime) {
    if (detections.length < this.config.minDetections) {
      return null;
    }

    // 初期震源を設定(最初に検知した観測点)
    const firstDetection = detections[0];
    let epicenter = {
      lat: Math.round(firstDetection.location[0] * 100) / 100,
      lon: Math.round(firstDetection.location[1] * 100) / 100,
      depth: this.config.initialDepth,
      originTime: firstDetection.detectionTime + this.config.initialTimeOffset * 1000
    };

    // Step 1: 大まかに震央を探索(深さ固定)
    epicenter = this.improveEpicenter(
      epicenter, detections, currentTime,
      this.config.gridStepLarge, this.config.gridStepLarge
    );

    // Step 2: 細かく震央を探索(深さ固定)
    epicenter = this.improveEpicenter(
      epicenter, detections, currentTime,
      this.config.gridStepSmall, this.config.gridStepSmall
    );

    // Step 3: 深さも含めて大まかに探索
    epicenter = this.improveEpicenter(
      epicenter, detections, currentTime,
      this.config.gridStepSmall, this.config.gridStepSmall,
      this.config.depthStepLarge
    );

    // Step 4: 深さも含めて細かく探索
    epicenter = this.improveEpicenter(
      epicenter, detections, currentTime,
      this.config.gridStepSmall, this.config.gridStepSmall,
      this.config.depthStepSmall
    );

    // 最終誤差レベルを計算
    const finalError = this.calculateErrorLevel(
      epicenter, detections, currentTime
    );

    this.currentEpicenter = {
      ...epicenter,
      errorLevel: finalError,
      numDetections: detections.length,
      timestamp: new Date(currentTime).toISOString()
    };

    return this.currentEpicenter;
  }
}

// === 統合検知システム ===
function detectEarthquakeWithEpicenter(
  intensityData,
  locations,
  historyMap = new Map(),
  activeEvents = new Map(),
  neighborCache = new Map(),
  detectionTimeMap = new Map(), // 検知時刻を記録
  epicenterEstimator = new EpicenterEstimator(),
  config = {
    threshold: 0.3,
    historyWindow: 10,
    minMagnitudeChange: 0.2,
    neighborRadius: 2,
    minNeighborCount: 1,
    neighborThreshold: 0.15,
    expirySeconds: {
      weaker: 15,
      weak: 30,
      medium: 60,
      strong: 120,
      stronger: 180
    },
    anomalyThreshold: 3.0,
    anomalyChangeLimit: 1.0,
    islandAnomalyThreshold: 5.0,
    enableEpicenterEstimation: true // 震源推定の有効化
  }
) {
  const now = Date.now();
  const detectedStations = new Set();
  const updatedEvents = new Map();
  const newDetections = [];

  // === 震度分類 ===
  const classifyIntensity = (intensity) => {
    if (intensity >= 5.0) return 'stronger';
    if (intensity >= 3.0) return 'strong';
    if (intensity >= 1.0) return 'medium';
    if (intensity >= -1.0) return 'weak';
    if (intensity >= -1.5) return 'weaker';
    return null;
  };

  // === 近隣観測点の取得 ===
  const getNeighbors = (index) => {
    const neighbors = [];
    const [lat, lon] = locations[index];

    for (let i = 0; i < locations.length; i++) {
      if (i === index) continue;

      const [nLat, nLon] = locations[i];
      const distance = Math.sqrt(
        Math.pow(lat - nLat, 2) + Math.pow(lon - nLon, 2)
      );

      if (distance <= config.neighborRadius) {
        neighbors.push(i);
      }
    }

    return neighbors;
  };

  // === 異常値判定 ===
  const isAnomaly = (stationId, intensity, history, isIsland = false) => {
    if (detectedStations.has(stationId)) return false;
    if (!history || history.length === 0) return false;

    const oldest = history[0]?.value ?? intensity;
    const change = Math.abs(intensity - oldest);
    const threshold = isIsland ? config.islandAnomalyThreshold : config.anomalyThreshold;

    return intensity >= threshold && change < config.anomalyChangeLimit;
  };

  // === Step 1: 履歴の更新と差分計算 ===
  const stationAnalysis = new Map();

  intensityData.forEach((intensity, index) => {
    const stationId = `station-${index}`;

    if (!historyMap.has(stationId)) {
      historyMap.set(stationId, []);
    }
    const history = historyMap.get(stationId);
    history.push({ time: now, value: intensity });

    while (history.length > 0 && now - history[0].time > config.historyWindow * 1000) {
      history.shift();
    }

    const oldest = history[0]?.value ?? intensity;
    const diff = intensity - oldest;

    stationAnalysis.set(stationId, {
      index,
      intensity,
      diff,
      history,
      location: locations[index]
    });
  });

  // === Step 2: 異常値の除外 ===
  const validStations = new Map();
  stationAnalysis.forEach((data, stationId) => {
    const isIsland = false;
    if (!isAnomaly(stationId, data.intensity, data.history, isIsland)) {
      validStations.set(stationId, data);
    }
  });

  // === Step 3: 揺れの検知 ===
  validStations.forEach((data, stationId) => {
    const { index, intensity, diff, location } = data;

    if (intensity < config.threshold || diff < config.minMagnitudeChange) {
      return;
    }

    if (!neighborCache.has(index)) {
      neighborCache.set(index, getNeighbors(index));
    }
    const neighbors = neighborCache.get(index);

    let activeNeighborCount = 0;
    for (const nIndex of neighbors) {
      const nStationId = `station-${nIndex}`;
      const nData = validStations.get(nStationId);

      if (nData && nData.diff >= config.neighborThreshold) {
        activeNeighborCount++;
      }
    }

    if (activeNeighborCount >= config.minNeighborCount) {
      detectedStations.add(stationId);
      
      // 検知時刻を記録
      if (!detectionTimeMap.has(stationId)) {
        detectionTimeMap.set(stationId, now);
      }

      // === Step 4: イベントの割り当て ===
      let assignedEventId = null;

      for (const nIndex of neighbors) {
        const nStationId = `station-${nIndex}`;

        for (const [eventId, event] of activeEvents) {
          if (event.stations.has(nStationId)) {
            if (!assignedEventId || event.createdAt < activeEvents.get(assignedEventId).createdAt) {
              assignedEventId = eventId;
            }
          }
        }
      }

      if (!assignedEventId) {
        assignedEventId = `event-${index}-${now}`;
      }

      if (!updatedEvents.has(assignedEventId)) {
        const existingEvent = activeEvents.get(assignedEventId);
        updatedEvents.set(assignedEventId, {
          id: assignedEventId,
          createdAt: existingEvent?.createdAt ?? now,
          lastUpdated: now,
          maxIntensity: existingEvent?.maxIntensity ?? intensity,
          classification: existingEvent?.classification ?? classifyIntensity(intensity),
          stations: new Set(existingEvent?.stations ?? []),
          expiresAt: new Map(existingEvent?.expiresAt ?? new Map()),
          epicenter: existingEvent?.epicenter ?? null // 震源情報を追加
        });
      }

      const event = updatedEvents.get(assignedEventId);
      event.lastUpdated = now;
      event.stations.add(stationId);

      if (intensity > event.maxIntensity) {
        const oldClassification = event.classification;
        event.maxIntensity = intensity;
        event.classification = classifyIntensity(intensity);

        if (event.classification !== oldClassification) {
          newDetections.push({
            eventId: assignedEventId,
            classification: event.classification,
            maxIntensity: intensity,
            upgraded: !!oldClassification
          });
        }
      }

      const classification = classifyIntensity(intensity);
      const expiryDuration = config.expirySeconds[classification] || config.expirySeconds.weak;
      event.expiresAt.set(stationId, now + expiryDuration * 1000);
    }
  });

  // === Step 5: 震源推定 ===
  if (config.enableEpicenterEstimation) {
    updatedEvents.forEach((event) => {
      // イベントに関連する検知情報を収集
      const eventDetections = [];
      event.stations.forEach((stationId) => {
        const data = stationAnalysis.get(stationId);
        const detectionTime = detectionTimeMap.get(stationId);
        if (data && detectionTime) {
          eventDetections.push({
            location: data.location,
            intensity: data.intensity,
            detectionTime: detectionTime
          });
        }
      });

      // 検知時刻順にソート
      eventDetections.sort((a, b) => a.detectionTime - b.detectionTime);

      // 震源を推定
      if (eventDetections.length >= 3) {
        const epicenter = epicenterEstimator.estimate(eventDetections, now);
        if (epicenter) {
          event.epicenter = epicenter;
        }
      }
    });
  }

  // === Step 6: イベントの期限管理 ===
  updatedEvents.forEach((event) => {
    event.expiresAt.forEach((expiryTime, stationId) => {
      if (now > expiryTime) {
        event.stations.delete(stationId);
        event.expiresAt.delete(stationId);
      }
    });

    if (event.stations.size === 0) {
      updatedEvents.delete(event.id);
    }
  });

  // === 戻り値の生成 ===
  const detectionResults = [];
  updatedEvents.forEach((event) => {
    event.stations.forEach((stationId) => {
      const data = stationAnalysis.get(stationId);
      if (data) {
        detectionResults.push({
          location: data.location,
          intensity: data.intensity,
          eventId: event.id,
          classification: event.classification,
          timestamp: new Date(now).toISOString(),
          epicenter: event.epicenter // 震源情報を含める
        });
      }
    });
  });

  return {
    detections: detectionResults,
    events: updatedEvents,
    newAlerts: newDetections,
    historyMap,
    neighborCache,
    detectionTimeMap,
    epicenterEstimator
  };
}