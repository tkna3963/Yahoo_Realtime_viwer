function now_time() {
    currentDateTime = new Date();
    return currentDateTime
}

function objectToJsonText(obj, space = 2) {
    try {
        // JSON.stringifyでオブジェクトをJSON形式に変換
        return JSON.stringify(obj, null, space);
    } catch (error) {
        // 変換エラーが発生した場合はエラーメッセージを返す
        return `Error converting object to JSON: ${error.message}`;
    }
}

var watch_id;

function Location_information_setup() {
    watch_id = navigator.geolocation.watchPosition(Location_information_update, function (e) { console.log(e.message); }, { "enableHighAccuracy": true, "timeout": 20000, "maximumAge": 2000 });
}

var Location_list = [];

function Location_information_update(position) {
    var latitude = position.coords.latitude;
    var longitude = position.coords.longitude;
    var accuracy = position.coords.accuracy;
    var altitudeAccuracy = position.coords.altitudeAccuracy;
    var heading = position.coords.heading;
    var speed = position.coords.speed;
    var result_data = [latitude, longitude, accuracy, altitudeAccuracy, heading, speed];
    Location_list.push(result_data);
    console.log(Location_list)
}

Location_information_setup()

function text_speak(test) {
    const utterance = new SpeechSynthesisUtterance(test)
    speechSynthesis.speak(utterance)
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
    const savedData = localStorage.getItem('siteListData');
    if (savedData) {
        return JSON.parse(savedData);
    } else {
        const xhr = new XMLHttpRequest();
        xhr.open('GET', 'https://weather-kyoshin.west.edge.storage-yahoo.jp/SiteList/sitelist.json', false); // 同期的にリクエストを行う
        xhr.send();
        if (xhr.status === 200) {
            const json_data = JSON.parse(xhr.responseText);
            localStorage.setItem('siteListData', JSON.stringify(json_data));
            return json_data;
        } else {
            console.error('Error loading settings:', xhr.status);
            return null;
        }
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

function calDistance(A, B) {
    var R = 6371;
    var dLat = deg2rad(B[0] - A[0]);
    var dLon = deg2rad(B[1] - A[1]);
    var a = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
        Math.cos(deg2rad(A[0])) * Math.cos(deg2rad(B[0])) *
        Math.sin(dLon / 2) * Math.sin(dLon / 2);
    var c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    var distance = R * c;
    return parseFloat(distance.toFixed(0));
}


function yahoo_Realtime_data() {
    now_time_date = now_time()
    currentDateTime = new Date();
    set_time = new Date(now_time_date.getTime() - 3 * 1000);
    url_time_date = Yahoo_Time_date_fromat(set_time)
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${url_time_date}.json`;
    //const apiUrl = "https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/20240417/20240417231516.json";
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();
    Request_status = `${xhr.status}${xhr.statusText}`
    const Yahoo_json_data = JSON.parse(xhr.responseText);
    format_get_date = formatYahooTimeDate(set_time)

    const strongEarthquake = Yahoo_json_data.realTimeData.intensity;

    if (Yahoo_json_data.hypoInfo === null) {
        yahoo_result_data_list = { "situation": "EEW_hasn't_been_issued", "get_date": format_get_date, "strongEarthquake": strongEarthquake, "url_time_date": url_time_date, "status": Request_status, "API_URL": apiUrl, "Telegram": Yahoo_json_data };
    } else {
        const reportId = Yahoo_json_data.hypoInfo.items[0].reportId;
        const reportNum = Yahoo_json_data.hypoInfo.items[0].reportNum;
        const reportTime = formatYahooTimeDate(Yahoo_json_data.hypoInfo.items[0].reportTime);
        const originTime = formatYahooTimeDate(Yahoo_json_data.hypoInfo.items[0].originTime);
        const set_timeDoriginTime = calculateTimeDifference(Yahoo_json_data.hypoInfo.items[0].originTime, set_time, "nomal");
        const set_timeDoriginTimeseconds = calculateTimeDifference(Yahoo_json_data.hypoInfo.items[0].originTime, set_time, "seconds");
        const regionName = Yahoo_json_data.hypoInfo.items[0].regionName;
        const calcIntensity = Yahoo_json_data.hypoInfo.items[0].calcintensity.replace("0", "");
        const magnitude = Yahoo_json_data.hypoInfo.items[0].magnitude;
        const depth = Yahoo_json_data.hypoInfo.items[0].depth.replace("km", "");
        const isFinal = Yahoo_json_data.hypoInfo.items[0].isFinal;
        const isTraining = Yahoo_json_data.hypoInfo.items[0].isTraining;
        const isCancel = Yahoo_json_data.hypoInfo.items[0].isCancel;
        const Wave_latitude = Yahoo_json_data.psWave.items[0].latitude.replace("N", "");
        const Wave_longitude = Yahoo_json_data.psWave.items[0].longitude.replace("E", "");
        let pRadius = parseFloat(Yahoo_json_data.psWave.items[0].pRadius).toFixed(0);
        let sRadius = parseFloat(Yahoo_json_data.psWave.items[0].sRadius).toFixed(0);
        yahoo_result_data_list = {
            "situation": "EEW_has_been_issued",
            "get_date": format_get_date,
            "strongEarthquake": strongEarthquake,
            "reportId": reportId,
            "reportNum": reportNum,
            "reportTime": reportTime,
            "originTime": originTime,
            "set_timeDoriginTime": set_timeDoriginTime,
            "set_timeDoriginTimeseconds": set_timeDoriginTimeseconds,
            "calcIntensity": calcIntensity,
            "regionName": regionName,
            "magnitude": magnitude,
            "depth": depth,
            "isFinal": isFinal,
            "isTraining": isTraining,
            "isCancel": isCancel,
            "Wave_latitude": Wave_latitude,
            "Wave_longitude": Wave_longitude,
            "pRadius": pRadius,
            "sRadius": sRadius,
            "status": Request_status,
            "url_time_date": url_time_date,
            "API_URL": apiUrl,
            "Telegram": Yahoo_json_data
        };
    }
    return yahoo_result_data_list
}