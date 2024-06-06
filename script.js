function now_time() {
    currentDateTime = new Date();
    return currentDateTime
}

function objectToText(obj) {
    let text = '';
    for (const [key, value] of Object.entries(obj)) {
        text += `"${key}":${value}, `;
    }
    return text.slice(0, -2);
}

function calculateTimeDifference(start_time_str, end_time_str) {
    var startTime = new Date(start_time_str);
    var endTime = new Date(end_time_str);
    var timeDifference = endTime - startTime;
    var seconds = Math.floor(timeDifference / 1000);
    var years = Math.floor(seconds / (3600 * 24 * 365.25));
    var months = Math.floor((seconds % (3600 * 24 * 365.25)) / (3600 * 24 * 30.4375));
    seconds %= (3600 * 24 * 30.4375);
    var days = Math.floor(seconds / (3600 * 24));
    seconds %= (3600 * 24);
    var hours = Math.floor(seconds / 3600);
    seconds %= 3600;
    var minutes = Math.floor(seconds / 60);
    seconds %= 60;

    if (years > 0) {
        return `${years}年${months}か月${days}日と${hours}時間${minutes}分${seconds}秒前`;
    } else if (months > 0) {
        return `${months}か月${days}日と${hours}時間${minutes}分${seconds}秒前`;
    } else if (days > 0) {
        return `${days}日と${hours}時間${minutes}分${seconds}秒前`;
    } else if (hours > 0) {
        return `${hours}時間${minutes}分${seconds}秒前`;
    } else if (minutes > 0) {
        return `${minutes}分${seconds}秒前`;
    } else {
        return `${seconds}秒前`;
    }
}

function Print_checker_for_debugging(data, name) {
    console.log(data, name)
}

function List_checker_for_debugging(listData, name) {
    var contentCounter = 1;
    var printText = '';
    for (var key in listData) {
        if (listData.hasOwnProperty(key)) {
            printText += `元リスト: ${name} ${contentCounter}番目のデータ: "${key}": ${listData[key]}\n\n`;
            contentCounter++;
        }
    }
    console.log(printText);
}

var num = 0;
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
}

Location_information_setup()

function text_speak(test) {
    const utterance = new SpeechSynthesisUtterance(test)
    speechSynthesis.speak(utterance)
}

function seismicIntensityConversion(char) {
    try {
        var unicodeValue = char.charCodeAt(0);
        var intensity = unicodeValue - 100;

        if (intensity === 1) {
            return -3;
        } else if (intensity === 2) {
            return -2.5;
        } else if (intensity === 3) {
            return -2;
        } else if (intensity === 4) {
            return -1.5;
        } else if (intensity === 5) {
            return -1;
        } else if (intensity === 6) {
            return -0.5;
        } else if (intensity === 7) {
            return 0;
        } else if (intensity === 8) {
            return 0.5;
        } else if (intensity === 9) {
            return 1;
        } else if (intensity === 10) {
            return 1.5;
        } else if (intensity === 11) {
            return 2;
        } else if (intensity === 12) {
            return 2.5;
        } else if (intensity === 13) {
            return 3;
        } else if (intensity === 14) {
            return 3.5;
        } else if (intensity === 15) {
            return 4;
        } else if (intensity === 16) {
            return 4.5;
        } else if (intensity === 17) {
            return 5;
        } else if (intensity === 18) {
            return 5.5;
        } else if (intensity === 19) {
            return 6;
        } else if (intensity === 20) {
            return 6.5;
        } else if (intensity === 21) {
            return 7;
        }

        return intensity;
    } catch (error) {
        console.error(error);
    }
}

function findMaxIntensityLocations(loc_list, seismicIntensityList) {
    if (loc_list.length !== seismicIntensityList.length) return console.error("Error: The lengths of loc_list and seismicIntensityList do not match.");
    let maxIntensity = -Infinity;
    let maxIntensityLocation = null;
    for (let i = 0; i < loc_list.length; i++) {
        let intensity = seismicIntensityList[i];
        let loc = loc_list[i];
        if (intensity > maxIntensity) {
            maxIntensity = intensity;
            maxIntensityLocation = loc;
        }
    }
    return maxIntensityLocation;
}

function findMinAndIndex(list) {
    var min = list[0];
    var minIndex = 0;
    for (var i = 1; i < list.length; i++) {
        if (list[i] < min) {
            min = list[i];
            minIndex = i;
        }
    }
    return { min: min, index: minIndex };
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

function findMaxIntensityLocations(loc_list, seismicIntensityList) {
    if (loc_list.length !== seismicIntensityList.length) return console.error("Error: The lengths of loc_list and seismicIntensityList do not match.");

    let maxIntensity = -Infinity;
    let maxIntensityLocation = null;

    for (let i = 0; i < loc_list.length; i++) {
        let intensity = seismicIntensityList[i];
        let loc = loc_list[i];

        if (intensity > maxIntensity) {
            maxIntensity = intensity;
            maxIntensityLocation = loc;
        }
    }

    return maxIntensityLocation;
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
    var timeMachineInput = document.getElementById("time_machine");
    var timeMachineValue = timeMachineInput.value;

    var set_time;
    if (timeMachineValue === "") {
        set_time_counter = 0;
        currentDateTime = new Date();
        set_time = new Date(now_time_date.getTime() - 3 * 1000);
    } else {
        set_time = new Date(timeMachineValue);
        set_time_counter += 1;
        set_time = new Date(set_time.getTime() + set_time_counter * 1000);
    }
    url_time_date = Yahoo_Time_date_fromat(set_time)
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${url_time_date}.json`;
    //const apiUrl = "https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/20240417/20240417231516.json";
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();
    Request_status = xhr.status
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
        const set_timeDoriginTime = calculateTimeDifference(Yahoo_json_data.hypoInfo.items[0].originTime, set_time);
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
        yahoo_result_data_list = { "situation": "EEW_has_been_issued", "get_date": format_get_date, "strongEarthquake": strongEarthquake, "reportId": reportId, "reportNum": reportNum, "reportTime": reportTime, "originTime": originTime, "set_timeDoriginTime": set_timeDoriginTime, "calcIntensity": calcIntensity, "regionName": regionName, "magnitude": magnitude, "depth": depth, "isFinal": isFinal, "isTraining": isTraining, "isCancel": isCancel, "Wave_latitude": Wave_latitude, "Wave_longitude": Wave_longitude, "pRadius": pRadius, "sRadius": sRadius, "status": Request_status, "url_time_date": url_time_date, "API_URL": apiUrl, "Telegram": Yahoo_json_data };
    }
    return yahoo_result_data_list
}

function surfaceGroundInformationProvisionAPI(現在地緯度, 現在地経度) {
    const xhr = new XMLHttpRequest();
    xhr.open('GET', `https://www.j-shis.bosai.go.jp/map/api/sstrct/V2/meshinfo.geojson?position=${現在地経度},${現在地緯度}&epsg=4301`, false); // 同期的にリクエストを行う
    xhr.send();
    if (xhr.status === 200) {
        const json_data = JSON.parse(xhr.responseText);
        return json_data;
    } else {
        console.error('Error loading settings:', xhr.status);
        return null;
    }
}

function deg2rad(deg) {
    return deg * (Math.PI / 180);
}

function distanceAttenuationType(necessitiesData) {
    const magJMA = necessitiesData[0];
    const depth = necessitiesData[1];
    const epicenterLocaltion = necessitiesData[2];
    const magW = magJMA - 0.171;
    const long = 10 ** (0.5 * magW - 1.85) / 2;
    const pointLocaltion = necessitiesData[3];
    const epicenterDistance = calDistance(epicenterLocaltion, pointLocaltion);
    const hypocenterDistance = (depth ** 2 + epicenterDistance ** 2) ** 0.5 - long;
    const x = Math.max(hypocenterDistance, 3);
    const gpv600 = 10 ** (
        0.58 * magW +
        0.0038 * depth - 1.29 -
        Math.log10(x + 0.0028 * (10 ** (0.5 * magW))) -
        0.002 * x
    );
    const arv = necessitiesData[4];
    const pgv400 = gpv600 * 1.31;
    const pgv = pgv400 * arv;
    const intensity = parseFloat((2.68 + 1.72 * Math.log10(pgv)).toFixed(2));
    var attenuationResultData = [intensity, epicenterDistance, parseFloat(pgv.toFixed(3)), arv];
    return attenuationResultData;
}

function parseStringArrayToNumberArray(stringArray) {
    return stringArray.map(str => parseFloat(str.replace(/"/g, '')));
}

let Device_info_list = {};
function Information_distribution_board() {
    const Earthquake_info_list = yahoo_Realtime_data()
    const Yahoo_locate_list = Yahoo_locate_data()["items"]
    Earthquake_info_list.Yahoo_locate_list = Yahoo_locate_list;
    Earthquake_info_list.seismicIntensityList = Earthquake_info_list["strongEarthquake"].split('').map(char => seismicIntensityConversion(char));
    Earthquake_info_list.MAXseismicIntensity = Math.max(...Earthquake_info_list["strongEarthquake"].split('').map(char => seismicIntensityConversion(char)));
    Earthquake_info_list.seismicIntensityListCount = Earthquake_info_list["strongEarthquake"].length
    Earthquake_info_list.Earthquake_intensity_greater_than_zero_list = Earthquake_info_list["seismicIntensityList"].filter(value => value >= 0);
    Earthquake_info_list.Earthquake_intensity_greater_than_zero_Count = Earthquake_info_list["seismicIntensityList"].filter(value => value >= 0).length;
    Earthquake_info_list.MaxIntensityLocation = findMaxIntensityLocations(Yahoo_locate_list, Earthquake_info_list["seismicIntensityList"])
    if (Location_list[Location_list.length - 1] !== undefined) {
        const Location_data = Location_list[Location_list.length - 1];
        Device_info_list.Current_location_latitude = Location_data[0]
        Device_info_list.Current_location_longitude = Location_data[1]
    } else {
        Device_info_list.Current_location_latitude = 35.658;
        Device_info_list.Current_location_longitude = 139.4413;
    }
    const yahoo_calDistance_list = [];
    for (var yahoo_locate_item of Yahoo_locate_list) {
        yahoo_calDistance_list.push(calDistance([Device_info_list["Current_location_latitude"], Device_info_list["Current_location_longitude"]], yahoo_locate_item));
    }
    Earthquake_info_list.nearest_seismograph = Yahoo_locate_list[findMinAndIndex(yahoo_calDistance_list)["index"]];
    Earthquake_info_list.nearest_earthquake = Earthquake_info_list["seismicIntensityList"][findMinAndIndex(yahoo_calDistance_list)["index"]];
    var surfaceGroundInformationProvision_data = surfaceGroundInformationProvisionAPI(Device_info_list["Current_location_latitude"], Device_info_list["Current_location_longitude"]);
    var necessities = [Earthquake_info_list["magnitude"], Earthquake_info_list["depth"], [Earthquake_info_list["Wave_latitude"], Earthquake_info_list["Wave_longitude"]], [Device_info_list["Current_location_latitude"], Device_info_list["Current_location_longitude"]], surfaceGroundInformationProvision_data.features[0].properties.ARV];
    var Distance_attenuation_type_data_list = distanceAttenuationType(necessities);
    Earthquake_info_list.distanceAttenuationType_calcIntensity = Distance_attenuation_type_data_list[0]
    Earthquake_info_list.epicenterDistance = Distance_attenuation_type_data_list[1]
    Earthquake_info_list.PGV = Distance_attenuation_type_data_list[2]
    Earthquake_info_list.ARV = Distance_attenuation_type_data_list[3]
    const result_data = { "Device_info_list": Device_info_list, "Earthquake_info_list": Earthquake_info_list };
    List_checker_for_debugging(result_data["Earthquake_info_list"])
    return result_data
}

function Message_Conversion(Telegram_List_data) {
    let text = '';
    const info = Telegram_List_data["Earthquake_info_list"];

    if (info["situation"] === "EEW_hasn't_been_issued") {
        text = `
    情報取得時刻:${info["get_date"]}
    リアルタイム最大強震震度:${info["MAXseismicIntensity"]}(${info["MaxIntensityLocation"]})
    現在の震度0以上の観測点数:${info["Earthquake_intensity_greater_than_zero_Count"]} (観測点全体の約${((info["Earthquake_intensity_greater_than_zero_Count"] / info["seismicIntensityListCount"]) * 100).toFixed(0)}% ${info["Earthquake_intensity_greater_than_zero_Count"]}/${info["seismicIntensityListCount"]})
    現在地付近の観測所の震度:${info["nearest_earthquake"]}(${info["nearest_seismograph"]})
 `;
    } else {
        text = `
    情報取得時刻:${info["get_date"]}
    リアルタイム最大強震震度:${info["MAXseismicIntensity"]}(${info["MaxIntensityLocation"]})                                                                                                                                                                                                                                                                  
    現在の震度0以上の観測点数:${info["Earthquake_intensity_greater_than_zero_Count"]} (観測点全体の約${((info["Earthquake_intensity_greater_than_zero_Count"] / info["seismicIntensityListCount"]) * 100).toFixed(0)}% ${info["Earthquake_intensity_greater_than_zero_Count"]}/${info["seismicIntensityListCount"]})
    現在地付近の観測所の震度:${info["nearest_earthquake"]}(${info["nearest_seismograph"]})
    
    緊急地震速報要素
    EventID:${info["reportId"]}
    通番:${info["reportNum"]}${info["isFinal"] === "true" ? "(最終報)" : ""}
    地震発生時刻:${info["originTime"]}(${info["set_timeDoriginTime"]})
    情報更新時刻:${info["reportTime"]}
    評価震央地点:${info["regionName"]} 現在地から約${info["epicenterDistance"]}km
    推定震度:${info["calcIntensity"]} 現在地予想震度:約${info["distanceAttenuationType_calcIntensity"].toFixed(0)}(${info["distanceAttenuationType_calcIntensity"]})
    推定規模（マグニチュード):${info["magnitude"]}
    推定深さ:${info["depth"]}km
    予想初期微動継続時間:約${((info["epicenterDistance"]/8).toFixed(0))}秒

    `;
    }
    return text;
}

function openImage(src) {
    window.open(src, "_blank", "width=600, height=400");
}
function openSubWindow() {
    window.open("Monitor.html", "_blank", "width=342px, height=440px");
}

function Kyoushin_moniter_url() {
    now_time_date = now_time()
    set_time = new Date(now_time_date.getTime() - 3 * 1000);
    url_time_date = Yahoo_Time_date_fromat(set_time)
    const Kyoushin_S_Url = `https://smi.lmoniexp.bosai.go.jp/data/map_img/RealTimeImg/jma_b/${url_time_date}.jma_b.gif`;
    const Tyoushuki_Url = `https://www.lmoni.bosai.go.jp/monitor/data/data/map_img/RealTimeImg/abrspmx_s/${url_time_date}.abrspmx_s.gif`;
    const PGA_Url = `https://smi.lmoniexp.bosai.go.jp/data/map_img/RealTimeImg/acmap_s/${url_time_date}.acmap_s.gif`
    return {
        "Kyoushin_S_Url": Kyoushin_S_Url,
        "Tyoushuki_Url": Tyoushuki_Url,
        "PGA_Url": PGA_Url
    };
}

function Update() {
    const urls = Kyoushin_moniter_url();
    document.getElementById('imgKyoushinS').src = urls.Kyoushin_S_Url;
    document.getElementById('imgTyoushuki').src = urls.Tyoushuki_Url;
    document.getElementById('imgPGA').src = urls.PGA_Url;
    console.log(urls);
}

function updateImages() {
    const sliderValue = document.getElementById('imageSlider').value;
    const images = document.querySelectorAll('.image-container');
    images.forEach((image, index) => {
        if (index < sliderValue) {
            image.style.display = 'block';
        } else {
            image.style.display = 'none';
        }
    });
}

document.getElementById('imageSlider').addEventListener('input', updateImages);

function MAP(datas) {
    var map = L.map('map').setView([datas["Device_info_list"]["Current_location_latitude"], datas["Device_info_list"]["Current_location_longitude"]], 8);
    L.tileLayer('https://api.maptiler.com/maps/jp-mierune-dark/256/{z}/{x}/{y}.png?key=IkxKT1ZrEb6IprnTOsui', {
        attribution: '&copy; <a href="https://www.maptiler.com/copyright/">MapTiler</a> contributors'
    }).addTo(map);
}

Update();
setInterval(Update, 1000);
updateImages();
