// === JMA2001走時表の実装（改良版） ===
class TravelTimeTable {
  constructor() {
    this.table = [];
    this.loaded = false;
    this.lookupCache = new Map(); // パフォーマンス向上用キャッシュ
  }

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

        this.table.push({ pTime, sTime, depth, distance });
      }
    }

    this.loaded = true;
    this.lookupCache.clear();
    console.log(`走時表データ読み込み完了: ${this.table.length}エントリ`);
  }

  setTableData(data) {
    this.table = data;
    this.loaded = true;
    this.lookupCache.clear();
  }

  linearInterpolate(x, x0, x1, y0, y1) {
    if (x1 === x0) return y0;
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }

  // バイリニア補間（改良版：キャッシュ付き）
  bilinearInterpolate(distance, depth, waveType = 'P') {
    // キャッシュキーを生成
    const cacheKey = `${distance.toFixed(2)}-${depth.toFixed(1)}-${waveType}`;
    if (this.lookupCache.has(cacheKey)) {
      return this.lookupCache.get(cacheKey);
    }

    if (!this.loaded || this.table.length === 0) {
      const avgVelocity = waveType === 'P' ? 7.0 : 4.0;
      return distance / avgVelocity + depth / 8.0;
    }

    const distances = [...new Set(this.table.map(e => e.distance))].sort((a, b) => a - b);
    const depths = [...new Set(this.table.map(e => e.depth))].sort((a, b) => a - b);

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

    const getPoint = (dist, dep) => {
      return this.table.find(e => e.distance === dist && e.depth === dep);
    };

    const p00 = getPoint(d0, z0);
    const p01 = getPoint(d0, z1);
    const p10 = getPoint(d1, z0);
    const p11 = getPoint(d1, z1);

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

      const result = waveType === 'P' ? nearest.pTime : nearest.sTime;
      this.lookupCache.set(cacheKey, result);
      return result;
    }

    const timeKey = waveType === 'P' ? 'pTime' : 'sTime';
    const t0 = this.linearInterpolate(distance, d0, d1, p00[timeKey], p10[timeKey]);
    const t1 = this.linearInterpolate(distance, d0, d1, p01[timeKey], p11[timeKey]);
    const result = this.linearInterpolate(depth, z0, z1, t0, t1);

    this.lookupCache.set(cacheKey, result);
    return result;
  }

  getPWaveTravelTime(epicentralDistance, depth) {
    return this.bilinearInterpolate(epicentralDistance, depth, 'P');
  }

  getSWaveTravelTime(epicentralDistance, depth) {
    return this.bilinearInterpolate(epicentralDistance, depth, 'S');
  }
}

const globalTravelTimeTable = new TravelTimeTable();

// 走時表データの拡充初期化
(function initializeTravelTimeTable() {
  const data = [];
  const depths = [0, 10, 30, 50, 100, 200, 300, 500, 700];
  const distances = [];
  
  // 距離グリッドを生成（近距離は密、遠距離は疎）
  for (let d = 0; d <= 50; d += 2) distances.push(d);
  for (let d = 55; d <= 100; d += 5) distances.push(d);
  for (let d = 110; d <= 300; d += 10) distances.push(d);
  for (let d = 320; d <= 600; d += 20) distances.push(d);
  for (let d = 650; d <= 1500; d += 50) distances.push(d);
  
  depths.forEach(depth => {
    distances.forEach(distance => {
      // P波速度モデル（簡易的な層構造を考慮）
      let pVelocity, sVelocity;
      
      if (depth < 30) {
        // 地殻
        pVelocity = 6.0 + depth * 0.01;
        sVelocity = pVelocity / 1.73;
      } else if (depth < 100) {
        // 上部マントル
        pVelocity = 7.8 + (depth - 30) * 0.005;
        sVelocity = pVelocity / 1.75;
      } else {
        // 下部マントル
        pVelocity = 8.1 + (depth - 100) * 0.002;
        sVelocity = pVelocity / 1.78;
      }
      
      // 震源距離を考慮した走時計算
      const hypocentralDist = Math.sqrt(distance * distance + depth * depth);
      let pTime, sTime;
      
      if (distance < 100) {
        // 直達波
        pTime = hypocentralDist / pVelocity;
        sTime = hypocentralDist / sVelocity;
      } else {
        // 屈折波を考慮（簡易モデル）
        const directTime = hypocentralDist / pVelocity;
        const refractedTime = depth / 6.5 + (distance - depth) / 8.0 + depth / 6.5;
        pTime = Math.min(directTime, refractedTime);
        sTime = pTime * (pVelocity / sVelocity);
      }
      
      data.push({ pTime, sTime, depth, distance });
    });
  });
  
  globalTravelTimeTable.setTableData(data);
})();

// === 震源推定クラス（大幅改良版） ===
class EpicenterEstimator {
  constructor(config = {}) {
    this.config = {
      initialDepth: 10,
      initialTimeOffset: -2,
      gridStepLarge: 0.5,
      gridStepSmall: 0.1,
      depthStepLarge: 50,
      depthStepSmall: 10,
      nearStationWeight: 1.0,
      nearStationRadius: 50,
      minDetections: 3,
      undetectedPenalty: 1.0,
      undetectedCheckTime: 3,
      useSWave: true, // S波を使用するかどうか
      multiStartOptimization: true, // 複数初期値から最適化
      numInitialPoints: 5, // 初期値の数
      convergenceThreshold: 0.001, // 収束判定閾値
      maxIterations: 50, // 最大反復回数
      ...config
    };
    
    this.currentEpicenter = null;
    this.estimationHistory = []; // デバッグ用履歴
    this.distanceCache = new Map(); // 距離計算キャッシュ
  }

  // 震央距離を計算（キャッシュ付き）
  calculateEpicentralDistance(lat1, lon1, lat2, lon2) {
    const key = `${lat1.toFixed(4)}-${lon1.toFixed(4)}-${lat2.toFixed(4)}-${lon2.toFixed(4)}`;
    if (this.distanceCache.has(key)) {
      return this.distanceCache.get(key);
    }

    const R = 6371;
    const dLat = (lat2 - lat1) * Math.PI / 180;
    const dLon = (lon2 - lon1) * Math.PI / 180;
    const a = Math.sin(dLat/2) * Math.sin(dLat/2) +
              Math.cos(lat1 * Math.PI / 180) * Math.cos(lat2 * Math.PI / 180) *
              Math.sin(dLon/2) * Math.sin(dLon/2);
    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
    const distance = R * c;

    this.distanceCache.set(key, distance);
    return distance;
  }

  // 方位角を計算
  calculateAzimuth(lat1, lon1, lat2, lon2) {
    const dLon = (lon2 - lon1) * Math.PI / 180;
    const y = Math.sin(dLon) * Math.cos(lat2 * Math.PI / 180);
    const x = Math.cos(lat1 * Math.PI / 180) * Math.sin(lat2 * Math.PI / 180) -
              Math.sin(lat1 * Math.PI / 180) * Math.cos(lat2 * Math.PI / 180) * Math.cos(dLon);
    return Math.atan2(y, x);
  }

  // 方位角カバレッジを計算（観測点が全方位に分布しているか）
  calculateAzimuthCoverage(epicenter, detections) {
    if (detections.length < 4) return 0.5; // 少ない場合は低評価

    const azimuths = detections.map(d => 
      this.calculateAzimuth(epicenter.lat, epicenter.lon, d.location[0], d.location[1])
    ).sort((a, b) => a - b);

    // 最大方位角ギャップを計算
    let maxGap = 0;
    for (let i = 0; i < azimuths.length; i++) {
      const gap = i < azimuths.length - 1 
        ? azimuths[i + 1] - azimuths[i]
        : 2 * Math.PI - azimuths[azimuths.length - 1] + azimuths[0];
      maxGap = Math.max(maxGap, gap);
    }

    // カバレッジスコア（ギャップが小さいほど高スコア）
    return 1.0 - (maxGap / (2 * Math.PI));
  }

  // 観測点重み付け（改良版）
  calculateStationWeight(epicenter, detection, allDetections) {
    const distance = this.calculateEpicentralDistance(
      epicenter.lat, epicenter.lon,
      detection.location[0], detection.location[1]
    );

    // 距離による重み（指数減衰）
    const distWeight = Math.exp(-distance / 200);

    // 方位角カバレッジ考慮
    const azimuthCoverage = this.calculateAzimuthCoverage(epicenter, allDetections);

    // 近傍観測点へのボーナス
    const nearBonus = distance < this.config.nearStationRadius ? 1.5 : 1.0;

    return distWeight * (0.5 + 0.5 * azimuthCoverage) * nearBonus;
  }

  // 誤差レベルを計算（P波+S波対応、改良版）
  calculateErrorLevel(epicenter, detections, currentTime) {
    const { lat, lon, depth, originTime } = epicenter;
    const residuals = [];
    const weights = [];
    let totalError = 0;

    detections.forEach(detection => {
      const distance = this.calculateEpicentralDistance(
        lat, lon, detection.location[0], detection.location[1]
      );

      // P波の残差
      const pTravelTime = globalTravelTimeTable.getPWaveTravelTime(distance, depth);
      const pCalculatedOriginTime = detection.detectionTime - pTravelTime * 1000;
      const pResidual = (pCalculatedOriginTime - originTime) / 1000; // 秒単位

      const weight = this.calculateStationWeight(epicenter, detection, detections);
      weights.push(weight);
      residuals.push(pResidual);

      totalError += (pResidual * pResidual) * weight;

      // S波も使用する場合（S波検知時刻があれば）
      if (this.config.useSWave && detection.sWaveTime) {
        const sTravelTime = globalTravelTimeTable.getSWaveTravelTime(distance, depth);
        const sCalculatedOriginTime = detection.sWaveTime - sTravelTime * 1000;
        const sResidual = (sCalculatedOriginTime - originTime) / 1000;

        // S波の重みは少し低めに設定
        const sWeight = weight * 0.7;
        totalError += (sResidual * sResidual) * sWeight;
      }
    });

    // 正規化
    const weightSum = weights.reduce((a, b) => a + b, 0);
    totalError = totalError / weightSum;

    // 深さに対する物理的制約（異常な深さにペナルティ）
    if (depth < 0) totalError += 1000;
    if (depth > 700) totalError += (depth - 700) * 0.1;

    return totalError;
  }

  // 震源候補を生成（改良版：対角方向も追加）
  generateCandidates(center, latStep, lonStep, depthStep = null) {
    const candidates = [];

    // 基本4方向
    candidates.push(
      { lat: center.lat + latStep, lon: center.lon, depth: center.depth, originTime: center.originTime },
      { lat: center.lat - latStep, lon: center.lon, depth: center.depth, originTime: center.originTime },
      { lat: center.lat, lon: center.lon + lonStep, depth: center.depth, originTime: center.originTime },
      { lat: center.lat, lon: center.lon - lonStep, depth: center.depth, originTime: center.originTime }
    );

    // 対角方向（より効率的な探索）
    candidates.push(
      { lat: center.lat + latStep, lon: center.lon + lonStep, depth: center.depth, originTime: center.originTime },
      { lat: center.lat + latStep, lon: center.lon - lonStep, depth: center.depth, originTime: center.originTime },
      { lat: center.lat - latStep, lon: center.lon + lonStep, depth: center.depth, originTime: center.originTime },
      { lat: center.lat - latStep, lon: center.lon - lonStep, depth: center.depth, originTime: center.originTime }
    );

    if (depthStep !== null) {
      candidates.push(
        { lat: center.lat, lon: center.lon, depth: center.depth + depthStep, originTime: center.originTime },
        { lat: center.lat, lon: center.lon, depth: Math.max(0, center.depth - depthStep), originTime: center.originTime }
      );
    }

    // 発生時刻の微調整
    const timeStep = 0.1; // 0.1秒単位
    candidates.push(
      { ...center, originTime: center.originTime + timeStep * 1000 },
      { ...center, originTime: center.originTime - timeStep * 1000 }
    );

    return candidates;
  }

  // 震源を反復的に改善（収束判定追加）
  improveEpicenter(initialEpicenter, detections, currentTime, latStep, lonStep, depthStep = null) {
    let currentEpicenter = { ...initialEpicenter };
    let improved = true;
    let iteration = 0;
    let previousError = this.calculateErrorLevel(currentEpicenter, detections, currentTime);

    while (improved && iteration < this.config.maxIterations) {
      improved = false;
      const currentError = this.calculateErrorLevel(currentEpicenter, detections, currentTime);

      const candidates = this.generateCandidates(currentEpicenter, latStep, lonStep, depthStep);

      for (const candidate of candidates) {
        const candidateError = this.calculateErrorLevel(candidate, detections, currentTime);

        if (candidateError < currentError) {
          currentEpicenter = { ...candidate };
          improved = true;

          // 収束判定
          if (Math.abs(candidateError - previousError) < this.config.convergenceThreshold) {
            improved = false;
          }
          previousError = candidateError;
          break;
        }
      }

      iteration++;
    }

    return currentEpicenter;
  }

  // 複数初期値からの最適化
  multiStartOptimization(detections, currentTime) {
    const results = [];

    // 初期値の候補を生成
    const initialPoints = [];
    
    // 1. 最初の検知点
    initialPoints.push({
      lat: detections[0].location[0],
      lon: detections[0].location[1],
      depth: this.config.initialDepth,
      originTime: detections[0].detectionTime + this.config.initialTimeOffset * 1000
    });

    // 2. 重心点
    const centerLat = detections.reduce((sum, d) => sum + d.location[0], 0) / detections.length;
    const centerLon = detections.reduce((sum, d) => sum + d.location[1], 0) / detections.length;
    initialPoints.push({
      lat: centerLat,
      lon: centerLon,
      depth: this.config.initialDepth,
      originTime: detections[0].detectionTime + this.config.initialTimeOffset * 1000
    });

    // 3. 各象限の代表点（検知が多い場合）
    if (detections.length >= 8) {
      const quadrants = [
        detections.filter(d => d.location[0] >= centerLat && d.location[1] >= centerLon),
        detections.filter(d => d.location[0] >= centerLat && d.location[1] < centerLon),
        detections.filter(d => d.location[0] < centerLat && d.location[1] >= centerLon),
        detections.filter(d => d.location[0] < centerLat && d.location[1] < centerLon)
      ];

      quadrants.forEach(quad => {
        if (quad.length > 0) {
          const qLat = quad.reduce((sum, d) => sum + d.location[0], 0) / quad.length;
          const qLon = quad.reduce((sum, d) => sum + d.location[1], 0) / quad.length;
          initialPoints.push({
            lat: qLat,
            lon: qLon,
            depth: this.config.initialDepth,
            originTime: detections[0].detectionTime + this.config.initialTimeOffset * 1000
          });
        }
      });
    }

    // 各初期値から最適化を実行
    initialPoints.slice(0, this.config.numInitialPoints).forEach(initial => {
      let epicenter = { ...initial };

      // 段階的な最適化
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepLarge, this.config.gridStepLarge);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall, this.config.depthStepLarge);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall, this.config.depthStepSmall);

      const error = this.calculateErrorLevel(epicenter, detections, currentTime);
      results.push({ epicenter, error });
    });

    // 最良の結果を返す
    results.sort((a, b) => a.error - b.error);
    return results[0].epicenter;
  }

  // 誤差統計を計算
  calculateErrorStatistics(epicenter, detections) {
    const residuals = [];
    
    detections.forEach(detection => {
      const distance = this.calculateEpicentralDistance(
        epicenter.lat, epicenter.lon,
        detection.location[0], detection.location[1]
      );

      const pTravelTime = globalTravelTimeTable.getPWaveTravelTime(distance, epicenter.depth);
      const pCalculatedOriginTime = detection.detectionTime - pTravelTime * 1000;
      const pResidual = Math.abs((pCalculatedOriginTime - epicenter.originTime) / 1000);
      residuals.push(pResidual);
    });

    residuals.sort((a, b) => a - b);
    const mean = residuals.reduce((a, b) => a + b, 0) / residuals.length;
    const variance = residuals.reduce((sum, r) => sum + Math.pow(r - mean, 2), 0) / residuals.length;
    const stdDev = Math.sqrt(variance);
    const median = residuals[Math.floor(residuals.length / 2)];

    // MAD (Median Absolute Deviation)
    const mad = residuals.map(r => Math.abs(r - median));
    mad.sort((a, b) => a - b);
    const madValue = mad[Math.floor(mad.length / 2)];

    return { mean, stdDev, median, mad: madValue, residuals };
  }

  // 外れ値検出
  detectOutliers(epicenter, detections) {
    const stats = this.calculateErrorStatistics(epicenter, detections);
    const outliers = [];

    detections.forEach((detection, index) => {
      const residual = stats.residuals[index];
      // MADベースの外れ値検出（3σ相当）
      if (Math.abs(residual - stats.median) > 3 * stats.mad) {
        outliers.push({ index, residual, detection });
      }
    });

    return { outliers, stats };
  }

  // 震源推定結果の妥当性検証
  validateEpicenter(epicenter, detections) {
    const warnings = [];
    let valid = true;

    // 深さの妥当性
    if (epicenter.depth < 0) {
      warnings.push('深さが負の値です');
      valid = false;
    }
    if (epicenter.depth > 700) {
      warnings.push('深さが異常に深い（>700km）');
    }

    // 観測点配置の妥当性
    const azimuthCoverage = this.calculateAzimuthCoverage(epicenter, detections);
    if (azimuthCoverage < 0.3) {
      warnings.push('観測点の方位角カバレッジが低い（偏った配置）');
    }

    // 誤差統計
    const { outliers, stats } = this.detectOutliers(epicenter, detections);
    if (outliers.length > detections.length * 0.3) {
      warnings.push(`外れ値が多い（${outliers.length}/${detections.length}）`);
    }

    if (stats.stdDev > 2.0) {
      warnings.push(`残差の標準偏差が大きい（${stats.stdDev.toFixed(2)}秒）`);
    }

    return { valid, warnings, stats, outliers, azimuthCoverage };
  }

  // メイン震源推定処理（大幅改良版）
  estimate(detections, currentTime) {
    if (detections.length < this.config.minDetections) {
      return null;
    }

    // 距離キャッシュをクリア（メモリ節約）
    if (this.distanceCache.size > 10000) {
      this.distanceCache.clear();
    }

    let epicenter;

    if (this.config.multiStartOptimization && detections.length >= 5) {
      // 複数初期値からの最適化
      epicenter = this.multiStartOptimization(detections, currentTime);
    } else {
      // 従来の単一初期値最適化
      const firstDetection = detections[0];
      epicenter = {
        lat: Math.round(firstDetection.location[0] * 100) / 100,
        lon: Math.round(firstDetection.location[1] * 100) / 100,
        depth: this.config.initialDepth,
        originTime: firstDetection.detectionTime + this.config.initialTimeOffset * 1000
      };

      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepLarge, this.config.gridStepLarge);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall, this.config.depthStepLarge);
      epicenter = this.improveEpicenter(epicenter, detections, currentTime,
        this.config.gridStepSmall, this.config.gridStepSmall, this.config.depthStepSmall);
    }

    // 誤差評価と妥当性検証
    const finalError = this.calculateErrorLevel(epicenter, detections, currentTime);
    const validation = this.validateEpicenter(epicenter, detections);

    this.currentEpicenter = {
      ...epicenter,
      errorLevel: finalError,
      numDetections: detections.length,
      timestamp: new Date(currentTime).toISOString(),
      validation: validation,
      confidence: this.calculateConfidence(validation)
    };

    // デバッグ履歴に追加
    this.estimationHistory.push({
      timestamp: currentTime,
      epicenter: this.currentEpicenter,
      numDetections: detections.length
    });

    // 履歴が長くなりすぎたら古いものを削除
    if (this.estimationHistory.length > 100) {
      this.estimationHistory.shift();
    }

    return this.currentEpicenter;
  }

  // 信頼度を計算（0-1のスコア）
  calculateConfidence(validation) {
    let confidence = 1.0;

    // 方位角カバレッジ
    confidence *= (0.5 + 0.5 * validation.azimuthCoverage);

    // 外れ値の割合
    const outlierRatio = validation.outliers.length / validation.stats.residuals.length;
    confidence *= (1.0 - outlierRatio);

    // 残差の大きさ
    confidence *= Math.exp(-validation.stats.stdDev / 2.0);

    return Math.max(0, Math.min(1, confidence));
  }
}

// === 統合検知システム（改良版） ===
function detectEarthquakeWithEpicenter(
  intensityData,
  locations,
  historyMap = new Map(),
  activeEvents = new Map(),
  neighborCache = new Map(),
  detectionTimeMap = new Map(),
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
    enableEpicenterEstimation: true,
    incrementalUpdate: true // 増分更新の有効化
  }
) {
  const now = Date.now();
  const detectedStations = new Set();
  const updatedEvents = new Map();
  const newDetections = [];

  const classifyIntensity = (intensity) => {
    if (intensity >= 5.0) return 'stronger';
    if (intensity >= 3.0) return 'strong';
    if (intensity >= 1.0) return 'medium';
    if (intensity >= -1.0) return 'weak';
    if (intensity >= -1.5) return 'weaker';
    return null;
  };

  const getNeighbors = (index) => {
    if (neighborCache.has(index)) {
      return neighborCache.get(index);
    }

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

    neighborCache.set(index, neighbors);
    return neighbors;
  };

  const isAnomaly = (stationId, intensity, history, isIsland = false) => {
    if (detectedStations.has(stationId)) return false;
    if (!history || history.length === 0) return false;

    const oldest = history[0]?.value ?? intensity;
    const change = Math.abs(intensity - oldest);
    const threshold = isIsland ? config.islandAnomalyThreshold : config.anomalyThreshold;

    return intensity >= threshold && change < config.anomalyChangeLimit;
  };

  // Step 1: 履歴の更新と差分計算
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

  // Step 2: 異常値の除外
  const validStations = new Map();
  stationAnalysis.forEach((data, stationId) => {
    const isIsland = false;
    if (!isAnomaly(stationId, data.intensity, data.history, isIsland)) {
      validStations.set(stationId, data);
    }
  });

  // Step 3: 揺れの検知
  validStations.forEach((data, stationId) => {
    const { index, intensity, diff, location } = data;

    if (intensity < config.threshold || diff < config.minMagnitudeChange) {
      return;
    }

    const neighbors = getNeighbors(index);

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
      
      if (!detectionTimeMap.has(stationId)) {
        detectionTimeMap.set(stationId, now);
      }

      // Step 4: イベントの割り当て
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
          epicenter: existingEvent?.epicenter ?? null,
          needsEpicenterUpdate: true // 震源更新フラグ
        });
      }

      const event = updatedEvents.get(assignedEventId);
      event.lastUpdated = now;
      event.stations.add(stationId);
      event.needsEpicenterUpdate = true;

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

  // Step 5: 震源推定（改良版：増分更新対応）
  if (config.enableEpicenterEstimation) {
    updatedEvents.forEach((event) => {
      // 増分更新が有効で、既に震源がある場合はスキップ可能
      if (config.incrementalUpdate && event.epicenter && !event.needsEpicenterUpdate) {
        return;
      }

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

      eventDetections.sort((a, b) => a.detectionTime - b.detectionTime);

      if (eventDetections.length >= 3) {
        const epicenter = epicenterEstimator.estimate(eventDetections, now);
        if (epicenter) {
          event.epicenter = epicenter;
          event.needsEpicenterUpdate = false;
        }
      }
    });
  }

  // Step 6: イベントの期限管理
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

  // 戻り値の生成
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
          epicenter: event.epicenter
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

// エクスポート（必要に応じて）
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    TravelTimeTable,
    EpicenterEstimator,
    detectEarthquakeWithEpicenter,
    globalTravelTimeTable
  };
}