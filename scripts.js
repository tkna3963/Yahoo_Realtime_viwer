let currentLocation = { latitude: null, longitude: null };

function updateLocation(position) {
    currentLocation.latitude = position.coords.latitude;
    currentLocation.longitude = position.coords.longitude;
    console.log(`現在地: 緯度 ${currentLocation.latitude}, 経度 ${currentLocation.longitude}`);
}

if (navigator.geolocation) {
    navigator.geolocation.watchPosition(updateLocation, error => {
        console.error('Geolocation error:', error);
    }, {
        enableHighAccuracy: true,
        maximumAge: 0,
        timeout: 1000
    });
} else {
    console.error("Geolocation is not supported by this browser.");
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

    // 月は手動で処理（2つ目の「年」が「月」の誤記に由来している前提p）
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

function getSeismicColor(intensity) {
    const colorMapping = {
        "-3": "#71a2cb", "-2.5": "#89afc8", "-2": "#90b3ca",
        "-1.5": "#97b7cc", "-1": "#97b7cc", "-0.5": "#97b7cc",
        "0": "#97b7cc", "0.5": "#90b3ca", "1": "#89afc8",
        "1.5": "#71a2cb", "2": "#5ea7ac", "2.5": "#38a477",
        "3": "#0fb02b", "3.5": "#f4e200", "4": "#fbc300",
        "4.5": "#ffaf00", "5": "#f90", "5.5": "#ff7e00",
        "6": "#ff6200", "6.5": "#fc4c02", "7": "#f53605"
    };

    return colorMapping[intensity] || "#000"; // 該当なしは黒
}

function getSeismicImage(intensity) {
    let result;
    const numIntensity = parseFloat(intensity); // intensityを数値に変換

    if (numIntensity < 0.5) {
        result = "Required_files/Shindo0.png"; // 震度0: 該当なし
    } else if (numIntensity < 1.5) {
        result = "Required_files/Shindo1.png"; // 震度1
    } else if (numIntensity < 2.5) {
        result = "Required_files/Shindo2.png"; // 震度2
    } else if (numIntensity < 3.5) {
        result = "Required_files/Shindo3.png"; // 震度3
    } else if (numIntensity < 4.5) {
        result = "Required_files/Shindo4.png"; // 震度4
    } else if (numIntensity < 5.0) {
        result = "Required_files/Shindo5-.png"; // 震度5弱
    } else if (numIntensity < 5.5) {
        result = "Required_files/Shindo5+.png"; // 震度5強
    } else if (numIntensity < 6.0) {
        result = "Required_files/Shindo6-.png"; // 震度6弱
    } else if (numIntensity < 6.5) {
        result = "Required_files/Shindo6+.png"; // 震度6強
    } else {
        result = "Required_files/Shindo7.png"; // 震度7
    }

    return result;
}

// 震度をラベルで返す関数（基準を未満に統一）
function getShindoLabel(shindo) {
    if (shindo < 0.5) return "0";
    if (shindo < 1.5) return "1";
    if (shindo < 2.5) return "2";
    if (shindo < 3.5) return "3";
    if (shindo < 4.5) return "4";
    if (shindo < 5.0) return "5弱";
    if (shindo < 5.5) return "5強";
    if (shindo < 6.0) return "6弱";
    if (shindo < 6.5) return "6強";
    return "7";
}

// 震度に対応する色を返す関数（基準を未満に統一）
function getColor(shindo) {
    let color;

    if (shindo < 0.5) {
        color = "#FFFFFF"; // 震度0: 白
    } else if (shindo < 1.5) {
        color = "#ADADAD"; // 震度1: グレー
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
    } else if (shindo <= 6.5) {
        color = "#5D1799"; // 震度7: 紫
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
        else if (value <= 6.5) shindo = "震度7";

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

function yahooRealtimeData() {
    var timeMachineInput = document.getElementById("time_machine");
    var timeMachineValue = timeMachineInput.value;

    if (timeMachineValue === "") {
        const currentDateTime = new Date();
        set_time = new Date(currentDateTime.getTime() - 3 * 1000);
    } else {
        set_time = new Date(timeMachineValue);
        set_time_counter += 1;
        set_time = new Date(set_time.getTime() + set_time_counter * 1000);
    }

    const url_time_date = Yahoo_Time_date_fromat(set_time); // Assuming this formats the date correctly.

    // Test URL for API (replace with `url_time_date` as needed)
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${url_time_date}.json`;

    // Make synchronous API request
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();

    const request_status = `${xhr.status}${xhr.statusText}`;
    const yahoo_json_data = JSON.parse(xhr.responseText);
    const format_get_date = formatYahooTimeDate(set_time); // Assuming this formats the date correctly.

    // Structure for no EEW
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

    // Structure for EEW issued
    const hypoInfo = yahoo_json_data.hypoInfo.items[0];
    const psWave = yahoo_json_data.psWave.items[0];

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
        pRadius: Math.round(psWave.pRadius) + 21,
        sRadius: Math.round(psWave.sRadius) + 12,
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
    const surfaceSpeed = maxSpeed400 * amplificationFactor;
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
        return cache.get(cacheKey);
    }
    loadKMData();
    let min_dist = Infinity;
    let closest_site = null;
    for (let km of KM_data) {
        let km_lat = parseFloat(km[3]);
        let km_lon = parseFloat(km[4]);
        if (isNaN(km_lat) || isNaN(km_lon)) continue;
        let dist = calculateDistance([latitude, longitude], [km_lat, km_lon]);
        if (dist < min_dist) {
            min_dist = dist;
            closest_site = km;
        }
    }
    if (closest_site) {
        let result = {
            latitude: latitude,
            longitude: longitude,
            closest_KM_site: closest_site[2],
            KM_lat: closest_site[3],
            KM_lon: closest_site[4],
            distance: min_dist
        };
        cache.set(cacheKey, result);
        return result;
    } else {
        console.log("近い地点が見つかりませんでした");
        return null;
    }
}

function datas_bord() {
    const results_datalist = {};
    const AreaSFClist = [];
    const Areanamelist = [];
    const Arvlist = [];
    const DBLlist = [];
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

        centers = loadJSON("Required_files/centerARVs.json")
        for (const center of centers.centers) {
            const AreaSFC = calculateDistanceAttenuation(
                YahooDatas.magnitude,
                YahooDatas.depth,
                [YahooDatas.Wave_latitude, YahooDatas.Wave_longitude],
                [center.latitude, center.longitude],
                Number(center.properties.arv)
            );
            AreaSFClist.push(AreaSFC.intensity);
            Areanamelist.push(center.properties.name);
            Arvlist.push(Number(center.properties.arv));
            DBLlist.push(AreaSFC.epicenterDistance);
        }
        results_datalist.AreaSFClist = AreaSFClist
        results_datalist.Areanamelist = Areanamelist
        results_datalist.Arvlist = Arvlist
        results_datalist.DBLlist = DBLlist
    }

    return results_datalist;
}

