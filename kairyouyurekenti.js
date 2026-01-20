// === オールインワン検知関数 ===
function newdetectEarthquakeEvents(
    intensityData,
    locations,
    historyMap = new Map(),
    activeEvents = new Map(),
    neighborCache = new Map(),
    config = {
        threshold: 0.5,
        historyWindow: 10,
        minMagnitudeChange: 0.5,
        neighborRadius: 2,
        minNeighborCount: 1,
        neighborThreshold: 0.3,
        expirySeconds: { weaker: 15, weak: 30, medium: 60, strong: 120, stronger: 180 },
        anomalyThreshold: 3.0,
        anomalyChangeLimit: 1.0,
        islandAnomalyThreshold: 5.0
    }
) {
    const now = Date.now();
    const detectedStations = new Set();
    const updatedEvents = new Map();
    const newDetections = [];

    // === 内部関数: 震度分類 ===
    const classifyIntensity = (intensity) => {
        if (intensity >= 5.0) return 'stronger';
        if (intensity >= 3.0) return 'strong';
        if (intensity >= 1.0) return 'medium';
        if (intensity >= -1.0) return 'weak';
        if (intensity >= -1.5) return 'weaker';
        return null;
    };

    // === 内部関数: 近隣観測点の取得 ===
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

    // === 内部関数: 異常値判定 ===
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

        // 履歴の更新
        if (!historyMap.has(stationId)) {
            historyMap.set(stationId, []);
        }
        const history = historyMap.get(stationId);
        history.push({ time: now, value: intensity });

        // 古い履歴の削除
        while (history.length > 0 && now - history[0].time > config.historyWindow * 1000) {
            history.shift();
        }

        // 差分計算
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
        const isIsland = false; // 実装時は離島判定ロジックを追加
        if (!isAnomaly(stationId, data.intensity, data.history, isIsland)) {
            validStations.set(stationId, data);
        }
    });

    // === Step 3: 揺れの検知 ===
    validStations.forEach((data, stationId) => {
        const { index, intensity, diff, location } = data;

        // 閾値チェック
        if (intensity < config.threshold || diff < config.minMagnitudeChange) {
            return;
        }

        // 近隣観測点の取得（キャッシュ利用）
        if (!neighborCache.has(index)) {
            neighborCache.set(index, getNeighbors(index));
        }
        const neighbors = neighborCache.get(index);

        // 近隣観測点の揺れをカウント
        let activeNeighborCount = 0;
        for (const nIndex of neighbors) {
            const nStationId = `station-${nIndex}`;
            const nData = validStations.get(nStationId);

            if (nData && nData.diff >= config.neighborThreshold) {
                activeNeighborCount++;
            }
        }

        // 揺れと判定
        if (activeNeighborCount >= config.minNeighborCount) {
            detectedStations.add(stationId);

            // === Step 4: イベントの割り当て ===
            let assignedEventId = null;

            // 近隣の既存イベントを探索
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

            // 新規イベント作成
            if (!assignedEventId) {
                assignedEventId = `event-${index}-${now}`;
            }

            // イベントの更新または作成
            if (!updatedEvents.has(assignedEventId)) {
                const existingEvent = activeEvents.get(assignedEventId);
                updatedEvents.set(assignedEventId, {
                    id: assignedEventId,
                    createdAt: existingEvent?.createdAt ?? now,
                    lastUpdated: now,
                    maxIntensity: existingEvent?.maxIntensity ?? intensity,
                    classification: existingEvent?.classification ?? classifyIntensity(intensity),
                    stations: new Set(existingEvent?.stations ?? []),
                    expiresAt: new Map(existingEvent?.expiresAt ?? new Map())
                });
            }

            const event = updatedEvents.get(assignedEventId);
            event.lastUpdated = now;
            event.stations.add(stationId);

            // 最大震度の更新
            if (intensity > event.maxIntensity) {
                const oldClassification = event.classification;
                event.maxIntensity = intensity;
                event.classification = classifyIntensity(intensity);

                // 震度分類が上昇した場合、通知が必要
                if (event.classification !== oldClassification) {
                    newDetections.push({
                        eventId: assignedEventId,
                        classification: event.classification,
                        maxIntensity: intensity,
                        upgraded: !!oldClassification
                    });
                }
            }

            // 観測点ごとの有効期限を設定
            const classification = classifyIntensity(intensity);
            const expiryDuration = config.expirySeconds[classification] || config.expirySeconds.weak;
            event.expiresAt.set(stationId, now + expiryDuration * 1000);
        }
    });

    // === Step 5: イベントの期限管理 ===
    updatedEvents.forEach((event) => {
        // 期限切れの観測点を削除
        event.expiresAt.forEach((expiryTime, stationId) => {
            if (now > expiryTime) {
                event.stations.delete(stationId);
                event.expiresAt.delete(stationId);
            }
        });

        // 観測点が0になったイベントは削除
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
                    timestamp: new Date(now).toISOString()
                });
            }
        });
    });

    return {
        detections: detectionResults,
        events: updatedEvents,
        newAlerts: newDetections,
        historyMap,
        neighborCache
    };
}

