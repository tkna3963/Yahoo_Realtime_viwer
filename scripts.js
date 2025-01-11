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
        timeout: 5000
    });
} else {
    console.error("Geolocation is not supported by this browser.");
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

function calDistance(A, B) {
    const R = 6371; // 地球の半径 (km)
    const dLat = deg2rad(B[0] - A[0]);
    const dLon = deg2rad(B[1] - A[1]);
    const a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
        Math.cos(deg2rad(A[0])) * Math.cos(deg2rad(B[0])) *
        Math.sin(dLon / 2) * Math.sin(dLon / 2);
    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    const distance = R * c;
    return parseFloat(distance.toFixed(2));
}

function yahoo_Realtime_data() {
    const now_time_date = now_time();
    const set_time = new Date(now_time_date.getTime() - 3 * 1000);
    const url_time_date = Yahoo_Time_date_fromat(set_time);
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${url_time_date}.json`;
    //const apiUrl = "https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/20240417/20240417231516.json";


    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();

    const Request_status = `${xhr.status}${xhr.statusText}`;
    const Yahoo_json_data = JSON.parse(xhr.responseText);
    const format_get_date = formatYahooTimeDate(set_time);

    if (Yahoo_json_data.hypoInfo === null) {
        return {
            "situation": "EEW_hasn't_been_issued",
            "get_date": format_get_date,
            "strongEarthquake": Yahoo_json_data.realTimeData.intensity,
            "url_time_date": url_time_date,
            "status": Request_status,
            "API_URL": apiUrl,
            "Telegram": Yahoo_json_data
        };
    } else {
        const info = Yahoo_json_data.hypoInfo.items[0];
        return {
            "situation": "EEW_has_been_issued",
            "get_date": format_get_date,
            "strongEarthquake": Yahoo_json_data.realTimeData.intensity,
            "reportId": info.reportId,
            "reportNum": info.reportNum,
            "reportTime": formatYahooTimeDate(info.reportTime),
            "originTime": formatYahooTimeDate(info.originTime),
            "regionName": info.regionName,
            "calcIntensity": info.calcintensity.replace("0", ""),
            "magnitude": info.magnitude,
            "depth": info.depth.replace("km", ""),
            "Wave_latitude": Yahoo_json_data.psWave.items[0].latitude.replace("N", ""),
            "Wave_longitude": Yahoo_json_data.psWave.items[0].longitude.replace("E", ""),
            "pRadius": Math.round(Yahoo_json_data.psWave.items[0].pRadius),
            "sRadius": Math.round(Yahoo_json_data.psWave.items[0].sRadius),
            "status": Request_status,
            "url_time_date": url_time_date,
            "API_URL": apiUrl,
            "Telegram": Yahoo_json_data
        };
    }
}

function datas_bord() {
    const results_datalist = {};
    const YahooDatas = yahoo_Realtime_data();

    results_datalist.yahoo_datas = YahooDatas;
    results_datalist.seismicIntensityList = YahooDatas["strongEarthquake"].split('').map(char => seismicIntensityConversion(char));
    results_datalist.MAXseismicIntensity = Math.max(...YahooDatas["strongEarthquake"].split('').map(char => seismicIntensityConversion(char)));
    results_datalist.latitude = currentLocation.latitude;
    results_datalist.longitude = currentLocation.longitude;

    return results_datalist;
}
