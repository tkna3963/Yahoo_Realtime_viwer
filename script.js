
let currentLocation = { latitude: 0, longitude: 0 };

// localStorage から保存された位置情報を読み込む関数
function loadLocationFromStorage() {
    const savedLocation = localStorage.getItem('savedLocation');
    if (savedLocation) {
        try {
            const parsed = JSON.parse(savedLocation);
            currentLocation.latitude = parsed.latitude;
            currentLocation.longitude = parsed.longitude;
            console.log(`ローカルストレージから読み込み: 緯度 ${currentLocation.latitude}, 経度 ${currentLocation.longitude}`);
        } catch (e) {
            console.error('ローカルストレージからの読み込みに失敗しました:', e);
        }
    } else {
        console.log('ローカルストレージに位置情報がありません。');
    }
}

function updateLocation(position) {
    currentLocation.latitude = position.coords.latitude;
    currentLocation.longitude = position.coords.longitude;
    console.log(`現在地: 緯度 ${currentLocation.latitude}, 経度 ${currentLocation.longitude}`);

    // 位置情報を localStorage に保存
    localStorage.setItem('savedLocation', JSON.stringify(currentLocation));
}

function handleLocationError(error) {
    console.error('Geolocation error:', error);

    // localStorage から保存された位置情報を読み込む
    loadLocationFromStorage();
}

if (navigator.geolocation) {
    navigator.geolocation.watchPosition(updateLocation, handleLocationError, {
        enableHighAccuracy: true,
        maximumAge: 0,
        timeout: 1000
    });
} else {
    console.error("このブラウザでは位置情報取得がサポートされていません。");
    loadLocationFromStorage(); // サポートされていない場合も試す
}

function playBeep(frequency = 440, duration = 500) {
    //console.log(frequency)
    const audioContext = new (window.AudioContext || window.webkitAudioContext)();
    const oscillator = audioContext.createOscillator();
    const gainNode = audioContext.createGain();

    oscillator.connect(gainNode);
    gainNode.connect(audioContext.destination);

    oscillator.type = 'sine'; // サイン波
    oscillator.frequency.setValueAtTime(frequency, audioContext.currentTime); // 周波数を設定
    gainNode.gain.setValueAtTime(1, audioContext.currentTime); // 音量を設定
    oscillator.start();

    setTimeout(() => {
        oscillator.stop();
    }, duration); // 指定した時間後に停止
}



function vibrate(pattern) {
    if ("vibrate" in navigator) {
        navigator.vibrate(pattern);
    } else {
        console.warn("このデバイスはバイブレーションをサポートしてないよ…");
    }
}

function notion(text) {
    // 通知の許可をリクエスト
    Notification.requestPermission().then((permission) => {
        if (permission === 'granted') {
            // 許可されている場合に通知を表示
            const notification = new Notification('通知タイトル', {
                body: text, // 通知の内容
            });
            // 通知がクリックされたときのイベントリスナー
            notification.onclick = () => {
                console.log('通知がクリックされました');
            };
        } else {
            console.log('通知が拒否されました');
        }
    });
}

function now_time() {
    const currentDateTime = new Date();
    return currentDateTime;
}

function objectToJsonText(obj, space = 2) {
    try {
        return JSON.stringify(obj, null, space);
    } catch (error) {
        return `Error converting object to JSON: ${error.message}`;
    }
}

const TemporaryDataRecordList = [];
function checkDuplicateData(Data) {
    TemporaryDataRecordList.push(Data);
    const isEqual = (a, b) => JSON.stringify(a) === JSON.stringify(b);
    const count = TemporaryDataRecordList.filter(
        data => isEqual(data, Data)
    ).length;
    if (count === 1) {
        return true;
    } else {
        if (count === 2) {
            const index = TemporaryDataRecordList.findIndex(
                data => isEqual(data, Data)
            );
            if (index !== -1) {
                TemporaryDataRecordList.splice(index, 1);
            }
        }
        return false;
    }
}

let previousTypes = {}; // reportIdごとの警報レベルを管理
function IsEEWType(data) {
    const reportId = data.reportId;
    let newType;
    if (data.calcIntensity === "6-" || data.calcIntensity === "6+" || data.calcIntensity === "7") {
        newType = "緊急地震速報(特別警報)";
    } else if (data.calcIntensity === "5-" || data.calcIntensity === "5+") {
        newType = "緊急地震速報(警報)";
    } else {
        newType = "緊急地震速報(予報)";
    }
    // reportIdごとに以前の警報レベルを管理
    if (!previousTypes[reportId]) {
        previousTypes[reportId] = newType;
        return newType;
    }
    const previousType = previousTypes[reportId];
    if (previousType === "緊急地震速報(警報)" && newType === "緊急地震速報(予報)") {
        return previousType;
    }
    // 更新
    previousTypes[reportId] = newType;
    return newType;
}

function Location_information_setup() {
    if (currentLocation.latitude && currentLocation.longitude) {
        console.log("位置情報が設定されました");
        return currentLocation;
    } else {
        console.error("位置情報が取得できていません");
        return null;
    }
}

function parseCustomDate(dateStr) {
    // 正規表現で日付と時刻を抽出（年の余分な「年」に対応）
    const match = dateStr.match(
        /(\d{4})年(?:\d{2})年(\d{2})日.*?(\d{2})時(\d{2})分(\d{2})秒/
    );

    if (!match) {
        console.error("日付のフォーマットが正しくありません:", dateStr);
        return NaN;
    }

    const [_, year, day, hour, minute, second] = match.map(Number);

    const month = parseInt(dateStr.match(/年(\d{2})年/)[1], 10);

    // Dateオブジェクトを作成
    return new Date(year, month - 1, day, hour, minute, second);
}


function calculateTimeDifference(reportTimeStr, eventTimeStr) {
    // 時間を文字列からDateオブジェクトに変換
    const reportTime = new Date(reportTimeStr);
    const eventTime = new Date(eventTimeStr);

    // 時差を計算（ミリ秒単位）
    const timeDifferenceMs = reportTime - eventTime;

    // ミリ秒を秒に変換
    const timeDifferenceSeconds = timeDifferenceMs / 1000;

    return timeDifferenceSeconds; // 秒単位の時差を返す
}

// 震央と地点の距離を計算するヘルパー関数（例: 緯度経度を使った距離計算）
function calculateDistance([lat1, lon1], [lat2, lon2]) {
    const R = 6371; // 地球の半径（km）
    const phi1 = lat1 * Math.PI / 180;
    const phi2 = lat2 * Math.PI / 180;
    const deltaPhi = (lat2 - lat1) * Math.PI / 180;
    const deltaLambda = (lon2 - lon1) * Math.PI / 180;

    const a = Math.sin(deltaPhi / 2) * Math.sin(deltaPhi / 2) +
        Math.cos(phi1) * Math.cos(phi2) *
        Math.sin(deltaLambda / 2) * Math.sin(deltaLambda / 2);

    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));

    return R * c; // キロメートル単位
}



function text_speak(test) {
    const utterance = new SpeechSynthesisUtterance(test);
    speechSynthesis.speak(utterance);
}

function seismicIntensityConversion(char) {
    try {
        const unicodeValue = char.charCodeAt(0);
        const intensity = unicodeValue - 100;

        const intensityMapping = {
            1: -3, 2: -2.5, 3: -2, 4: -1.5, 5: -1,
            6: -0.5, 7: 0, 8: 0.5, 9: 1, 10: 1.5,
            11: 2, 12: 2.5, 13: 3, 14: 3.5, 15: 4,
            16: 4.5, 17: 5, 18: 5.5, 19: 6, 20: 6.5,
            21: 7
        };

        return intensityMapping[intensity] ?? intensity;
    } catch (error) {
        console.error(error);
    }
}



// 震度をラベルで返す関数（基準を未満に統一）
function getShindoLabel(shindo) {
    if (shindo < -3) return "-3"; // 震度-3
    if (shindo < -2.5) return "-2.5"; // 震度-2.5
    if (shindo < -2) return "-2"; // 震度-2
    if (shindo < -1.5) return "-1.5"; // 震度-1.5
    if (shindo < -1.0) return "-1"; // 震度-1
    if (shindo < -0.5) return "-0.5"; // 震度-0.5
    if (shindo < 0.5) return "0"; // 震度0
    if (shindo < 1.5) return "1"; // 震度1
    if (shindo < 2.5) return "2"; // 震度2
    if (shindo < 3.5) return "3"; // 震度3
    if (shindo < 4.5) return "4"; // 震度4
    if (shindo < 5.0) return "5弱"; // 震度5弱
    if (shindo < 5.5) return "5強"; // 震度5強
    if (shindo < 6.0) return "6弱"; // 震度6弱
    if (shindo < 6.5) return "6強"; // 震度6強
    if (shindo <= 7.0) return "7"; // 震度7
    return "不明"; // それ以上
}

// 震度に対応する色を返す関数（基準を未満に統一）
function getSeismicColor(shindo) {
    let color;
    if (shindo < -3) {
        color = "#97b7cc"; // 震度-3:
    } else if (shindo < -2.5) {
        color = "#90b3ca"; // 震度-2.5:
    } else if (shindo < -2) {
        color = "#89afc8"; // 震度-2以下:
    } else if (shindo < -1.5) {
        color = "#71a2cb"; // 震度-1.5:
    } else if (shindo < -1.0) {
        color = "#5ea7ac"; // 震度-1:
    } else if (shindo < -0.5) {
        color = "#38a477"; // 震度-0.5:
    } else if (shindo < 0.5) {
        color = "#ADADAD"; // 震度0: グレー
    } else if (shindo < 1.5) {
        color = "#FFFFFF"; // 震度1: 白
    } else if (shindo < 2.5) {
        color = "#0066CC"; // 震度2: 青
    } else if (shindo < 3.5) {
        color = "#00A652"; // 震度3: 緑
    } else if (shindo < 4.5) {
        color = "#E6C300"; // 震度4: 黄色/金
    } else if (shindo < 5.0) {
        color = "#DC6B00"; // 震度5弱: オレンジ
    } else if (shindo < 5.5) {
        color = "#DC6B00"; // 震度5強: オレンジ (同じ色で表示)
    } else if (shindo < 6.0) {
        color = "#E33977"; // 震度6弱: ピンク/マゼンタ
    } else if (shindo < 6.5) {
        color = "#E33977"; // 震度6強: ピンク/マゼンタ (同じ色で表示)
    } else if (shindo <= 7.0) {
        color = "#5D1799"; // 震度7: 紫
    } else {
        color = "#1c0142"; // それ以上 (未知の震度)
    }
    return color;
}


// 震度データを受け取り、震度別に地域数をカウントする関数
function countSeismicIntensity(data) {
    const bins = {
        "震度7": 0, "震度6強": 0, "震度6弱": 0, "震度5強": 0,
        "震度5弱": 0, "震度4": 0, "震度3": 0,
        "震度2": 0, "震度1": 0
    };

    data.forEach(value => {
        let shindo;
        if (value < 0.5) shindo = "震度0";
        else if (value < 1.5) shindo = "震度1";
        else if (value < 2.5) shindo = "震度2";
        else if (value < 3.5) shindo = "震度3";
        else if (value < 4.5) shindo = "震度4";
        else if (value < 5.0) shindo = "震度5弱";
        else if (value < 5.5) shindo = "震度5強";
        else if (value < 6.0) shindo = "震度6弱";
        else if (value < 6.5) shindo = "震度6強";
        else if (6.5 <= value) shindo = "震度7";

        bins[shindo]++;
    });

    return Object.entries(bins)
        .filter(([_, count]) => count > 0) // 0件の震度を除外
        .map(([key, count]) => ` --${key}:${count}地域--\n\n`);
}


function Yahoo_locate_data() {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', 'https://weather-kyoshin.west.edge.storage-yahoo.jp/SiteList/sitelist.json', false); // 同期的リクエスト
    xhr.send();
    if (xhr.status === 200) {
        const json_data = JSON.parse(xhr.responseText);
        return json_data;
    } else {
        console.error('Error loading settings:', xhr.status);
        return null;
    }
}

function padZero(value) {
    return value < 10 ? '0' + value : value;
}

function Yahoo_Time_date_fromat(dateTime) {
    const year = dateTime.getFullYear();
    const month = padZero(dateTime.getMonth() + 1);
    const date = padZero(dateTime.getDate());
    const hours = padZero(dateTime.getHours());
    const minutes = padZero(dateTime.getMinutes());
    const seconds = padZero(dateTime.getSeconds());
    return `${year}${month}${date}/${year}${month}${date}${hours}${minutes}${seconds}`;
}

function formatYahooTimeDate(dateString, timeZone = 'Asia/Tokyo') {
    const options = { year: 'numeric', month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit', second: '2-digit', timeZone: timeZone };
    const date = new Date(dateString);
    const formatter = new Intl.DateTimeFormat('ja-JP', options);
    const formattedDate = formatter.format(date);
    const dayOfWeek = ['日', '月', '火', '水', '木', '金', '土'][date.getDay()];
    return formattedDate.replace(/\//g, '年').replace(' ', `日(${dayOfWeek}曜日) `).replace(':', '時').replace(':', '分') + '秒';
}

function deg2rad(deg) {
    return deg * (Math.PI / 180);
}

let set_time_counter = 0;
let set_time;

let past_pl = [];
let past_sl = [];

function yahooRealtimeData() {
    var timeMachineInput = document.getElementById("time_machine");
    var timeMachineValue = timeMachineInput.value;

    if (timeMachineValue === "") {
        const currentDateTime = new Date();
        set_time = new Date(currentDateTime.getTime() - 3 * 1000);
        document.getElementById("get_time").style.color = "white";
    } else {
        set_time = new Date(timeMachineValue);
        set_time_counter += 1;
        set_time = new Date(set_time.getTime() + set_time_counter * 1000);
        document.getElementById("get_time").style.color = "yellow";
    }

    const url_time_date = Yahoo_Time_date_fromat(set_time); // 日付フォーマット関数
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${url_time_date}.json`;

    // 同期リクエスト（同期の必要がないなら非同期に変更推奨）
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();

    const request_status = `${xhr.status}${xhr.statusText}`;
    const yahoo_json_data = JSON.parse(xhr.responseText);
    const format_get_date = formatYahooTimeDate(set_time); // 表示用日付フォーマット

    // 緊急地震速報が発表されていない場合
    if (yahoo_json_data.hypoInfo === null) {
        return {
            situation: "EEW_hasn't_been_issued",
            get_date: format_get_date,
            strongEarthquake: yahoo_json_data.realTimeData.intensity,
            url_time_date: url_time_date,
            status: request_status,
            API_URL: apiUrl,
            Telegram: yahoo_json_data
        };
    }

    // 緊急地震速報がある場合
    const hypoInfo = yahoo_json_data.hypoInfo.items[0];
    const psWave = yahoo_json_data.psWave.items[0];

    const pRadius = Math.round(psWave.pRadius);
    const sRadius = Math.round(psWave.sRadius);

    // 差分計算（前回値がある場合）
    let past_p_diff = null;
    let past_s_diff = null;

    if (past_pl.length > 0) {
        past_p_diff = pRadius - past_pl[past_pl.length - 1];
    }
    if (past_sl.length > 0) {
        past_s_diff = sRadius - past_sl[past_sl.length - 1];
    }

    // 履歴に現在値を追加
    past_pl.push(pRadius);
    past_sl.push(sRadius);

    return {
        situation: "EEW_has_been_issued",
        get_date: format_get_date,
        strongEarthquake: yahoo_json_data.realTimeData.intensity,
        reportId: hypoInfo.reportId,
        reportNum: hypoInfo.reportNum,
        isFinal: hypoInfo.isFinal,
        reportTime: formatYahooTimeDate(hypoInfo.reportTime),
        originTime: formatYahooTimeDate(hypoInfo.originTime),
        regionName: hypoInfo.regionName,
        calcIntensity: hypoInfo.calcintensity.replace("0", ""),
        magnitude: hypoInfo.magnitude,
        depth: hypoInfo.depth.replace("km", ""),
        Wave_latitude: psWave.latitude.replace("N", ""),
        Wave_longitude: psWave.longitude.replace("E", ""),
        pRadius: pRadius,
        sRadius: sRadius,
        past_p_diff: past_p_diff,
        past_s_diff: past_s_diff,
        status: request_status,
        url_time_date: url_time_date,
        API_URL: apiUrl,
        Telegram: yahoo_json_data
    };
}


// 地表情報提供APIを呼び出す関数
function surfaceGroundInformationProvisionAPI(現在地緯度, 現在地経度) {
    // ローカルストレージからデータを取得
    const savedData = localStorage.getItem('surfaceARV');
    // ローカルストレージに保存されたデータがあればそれを返す
    if (savedData) {
        return JSON.parse(savedData);
    } else {
        // XMLHttpRequestを使用して同期的にJSONファイルを読み込む
        const xhr = new XMLHttpRequest();
        xhr.open('GET', `https://www.j-shis.bosai.go.jp/map/api/sstrct/V2/meshinfo.geojson?position=${現在地経度},${現在地緯度}&epsg=4301`, false); // 同期的にリクエストを行う
        xhr.send();
        // レスポンスが成功（ステータスコード200）の場合のみ処理
        if (xhr.status === 200) {
            const json_data = JSON.parse(xhr.responseText);
            // ローカルストレージにデータを保存
            localStorage.setItem('surfaceARV', JSON.stringify(json_data));
            return json_data; // JSONを解析して返す
        } else {
            console.error('Error loading settings:', xhr.status);
            return null;
        }
    }
}
// 距離減衰タイプを計算する関数
function calculateDistanceAttenuation(magJMA, depth, epicenterLocation, pointLocation, amplificationFactor) {
    // マグニチュード (Mw) の計算
    const magW = magJMA - 0.171;
    // 断層長（半径）の計算
    const faultRadius = 10 ** (0.5 * magW - 1.85) / 2;
    // 震央からの距離の計算
    const epicenterDistance = calculateDistance(epicenterLocation, pointLocation);
    // 震源からの距離を計算（最小値3km）
    const hypocenterDistance = Math.max(Math.sqrt(depth ** 2 + epicenterDistance ** 2) - faultRadius, 3);
    // 工学基盤上の最大速度（Vs = 600 m/s）を計算
    const maxSpeed600 = 10 ** (
        0.58 * magW +
        0.0038 * depth - 1.29 -
        Math.log10(hypocenterDistance + 0.0028 * 10 ** (0.5 * magW)) -
        0.002 * hypocenterDistance
    );
    // 基盤の速度（Vs = 400 m/s）の変換
    const maxSpeed400 = maxSpeed600 * 1.31;
    // 増幅率を使って最終的な速度を計算
    const surfaceSpeed = amplificationFactor === 0 ? 0 : maxSpeed400 * (Number(amplificationFactor));
    // 震度を計算
    const intensity = parseFloat((2.68 + 1.72 * Math.log10(surfaceSpeed)).toFixed(2));
    // 結果を返す（震度、震央からの距離、最大速度、増幅率）
    return { intensity, epicenterDistance, surfaceSpeed: parseFloat(surfaceSpeed.toFixed(3)), amplificationFactor };
}

function loadJSON(filePath) {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', filePath, false);
    xhr.send();
    if (xhr.status === 200) {
        return JSON.parse(xhr.responseText);
    } else {
        console.error('Error loading JSON:', xhr.status, xhr.statusText);
        return null;
    }
}

function loadCSV(filePath) {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', filePath, false); // 同期リクエスト
    xhr.send();
    if (xhr.status === 200) {
        return xhr.responseText.split('\n').map(row => row.split(','));; // そのまま返す
    } else {
        console.error('Error loading CSV:', xhr.status, xhr.statusText);
        return null;
    }
}

const cache = new Map(); // キャッシュ用
let KM_data = null;
function loadKMData() {
    if (!KM_data) {
        KM_data = loadCSV("Required_files/sitepub_all_sj.csv");
    }
}

function YM_KM_C(latitude, longitude) {
    let cacheKey = `${latitude},${longitude}`;
    if (cache.has(cacheKey)) {
        return cache.get(cacheKey);  // キャッシュがある場合はそれを返す
    }
    loadKMData();
    let min_dist = Infinity;
    let closest_site = null;
    // 同じ場所が何度も選ばれないようにリストから除外
    let visitedSites = new Set();  // 訪れたサイトを追跡
    for (let km of KM_data) {
        let km_lat = parseFloat(km[4]);
        let km_lon = parseFloat(km[5]);
        if (isNaN(km_lat) || isNaN(km_lon)) continue;
        // すでに訪れたサイトならスキップ
        if (visitedSites.has(km[1])) {
            continue;
        }
        let dist = calculateDistance([latitude, longitude], [km_lat, km_lon]);
        if (dist < min_dist) {
            min_dist = dist;
            closest_site = km;
        }
        visitedSites.add(km[1]);  // 訪れたサイトを記録
    }
    if (closest_site) {
        let result = {
            id: closest_site[1],
            latitude: latitude,
            longitude: longitude,
            closest_KM_site: closest_site[2],
            KM_lat: closest_site[4],
            KM_lon: closest_site[5],
            distance: min_dist.toFixed(2)
        };
        cache.set(cacheKey, result);  // キャッシュに格納
        return result;
    } else {
        console.log("近い地点が見つかりませんでした");
        return null;
    }
}

let CentersData = loadJSON("Required_files\\SaibuncenterARVs.json")
function datas_bord() {
    const results_datalist = {};
    const AreaSFClist = [];
    const Areanamelist = [];
    const Arvlist = [];
    const DBLlist = [];
    const centerlalolist=[];
    const YahooDatas = yahooRealtimeData();
    results_datalist.yahoo_datas = YahooDatas;
    results_datalist.seismicIntensityList = YahooDatas["strongEarthquake"].split('').map(char => seismicIntensityConversion(char));
    results_datalist.MAXseismicIntensity = Math.max(...YahooDatas["strongEarthquake"].split('').map(char => seismicIntensityConversion(char)));
    results_datalist.latitude = currentLocation.latitude;
    results_datalist.longitude = currentLocation.longitude;
    if (currentLocation.latitude, currentLocation.longitude, YahooDatas.magnitude) {
        var SFGIPJ = surfaceGroundInformationProvisionAPI(currentLocation.latitude, currentLocation.longitude)
        var SFGIPARV = SFGIPJ.features[0].properties.ARV
        var SFC = calculateDistanceAttenuation(YahooDatas.magnitude, YahooDatas.depth, [YahooDatas.Wave_latitude, YahooDatas.Wave_longitude], [currentLocation.latitude, currentLocation.longitude], SFGIPARV)
        var SFCI = SFC.intensity
        var SFCD = Math.round(SFC.epicenterDistance)
        results_datalist.SFCI = SFCI;
        results_datalist.SFCD = SFCD;
        for (const centerP of CentersData.centers) {
            const AreaSFC = calculateDistanceAttenuation(
                YahooDatas.magnitude,
                YahooDatas.depth,
                [YahooDatas.Wave_latitude, YahooDatas.Wave_longitude],
                [centerP.latitude, centerP.longitude],
                Number(centerP.properties.arv)
            );
            AreaSFClist.push(AreaSFC.intensity);
            Areanamelist.push(centerP.properties.name);
            Arvlist.push(Number(centerP.properties.arv));
            DBLlist.push(AreaSFC.epicenterDistance);
            centerlalolist.push([centerP.latitude,centerP.longitude])
        }
        results_datalist.AreaSFClist = AreaSFClist
        results_datalist.Areanamelist = Areanamelist
        results_datalist.Arvlist = Arvlist
        results_datalist.DBLlist = DBLlist
        results_datalist.centerlalolist=centerlalolist
    }

    return results_datalist;
}


const threshold = 0.5; // 揺れとみなすリアルタイム震度
const historyWindow = 10; // 過去10秒間
const minMagnitudeChange = 0.5; // 震度の変化量での検知基準
const eventExpirySeconds = 15; // イベントを継続する秒数

function detectEvents(intensityData, locations, historyMap = new Map(), previousEvents = new Map()) {
    const now = Date.now();
    const events = [];
    const newEventMap = new Map();
    const eventAssignments = new Map(); // 観測点 -> イベント

    intensityData.forEach((intensity, index) => {
        const location = locations[index];
        const id = `station-${index}`;

        // === 観測点ごとの履歴を保持 ===
        if (!historyMap.has(id)) historyMap.set(id, []);
        const history = historyMap.get(id);
        history.push({ time: now, value: intensity });
        // 過去10秒より古い履歴を削除
        while (history.length > 0 && now - history[0].time > historyWindow * 1000) {
            history.shift();
        }

        // === 差分を計算 ===
        const oldest = history[0]?.value ?? intensity;
        const diff = intensity - oldest;

        // === 閾値を超えていたら処理開始 ===
        if (intensity >= threshold && diff >= minMagnitudeChange) {
            // 近隣インデックスを取得
            const neighbors = [index - 1, index + 1].filter(i => i >= 0 && i < intensityData.length);
            let mergedEventId = null;

            // 近隣の観測点が既存イベントに属しているか確認
            for (const ni of neighbors) {
                const neighborId = `station-${ni}`;
                const eventId = eventAssignments.get(neighborId);
                if (eventId && previousEvents.has(eventId)) {
                    mergedEventId = eventId;
                    break;
                }
            }

            // イベントIDを決定
            const eventId = mergedEventId || `event-${index}-${now}`;

            // 既存イベントがあるなら更新、ないなら新規作成
            if (!newEventMap.has(eventId)) {
                newEventMap.set(eventId, {
                    id: eventId,
                    createdAt: now,
                    lastUpdated: now,
                    maxIntensity: intensity,
                    stationIds: new Set([id]),
                });
            } else {
                const ev = newEventMap.get(eventId);
                ev.lastUpdated = now;
                ev.maxIntensity = Math.max(ev.maxIntensity, intensity);
                ev.stationIds.add(id);
            }

            // 観測点をイベントに紐づけ
            eventAssignments.set(id, eventId);

            // === 返却形式に合わせてイベントを追加 ===
            const currentTime = new Date(now).toISOString();
            events.push([location[0], location[1], currentTime]);
        }
    });

    return events;
}

