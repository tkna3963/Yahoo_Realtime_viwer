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

function padZero(value) {
    return value < 10 ? '0' + value : value;
}

function formatDateTimeForUrl(dateTime) {
    const year = dateTime.getFullYear();
    const month = padZero(dateTime.getMonth() + 1);
    const date = padZero(dateTime.getDate());
    const hours = padZero(dateTime.getHours());
    const minutes = padZero(dateTime.getMinutes());
    const seconds = padZero(dateTime.getSeconds());
    return `${year}${month}${date}/${year}${month}${date}${hours}${minutes}${seconds}`;
}

function yahooShingenn() {
    //const currentDateTime = new Date();
    const currentDateTime = new Date(2024, 0, 1, 16, 6, 19);
    const fiveSecondsAgo = new Date(currentDateTime.getTime() - 5 * 1000);
    const time_set = formatDateTimeForUrl(fiveSecondsAgo);
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${time_set}.json`;
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();

    if (xhr.readyState === 4 && xhr.status === 200) {
        const yahoo_data = JSON.parse(xhr.responseText);
        if (yahoo_data.hypoInfo === null) {
            const strongEarthquake = yahoo_data.realTimeData.intensity;
            const maxstrongEarthquake = Math.max(...strongEarthquake.split('').map(char => seismicIntensityConversion(char)));
            return [currentDateTime,maxstrongEarthquake, apiUrl, 0, "ホームポイント","緊急地震速報はでていません。", 0, 0,36.0,137.9, 0, 0];
        } else {
            const reportNum = yahoo_data.hypoInfo.items[0].reportNum;
            const isFinal = yahoo_data.hypoInfo.items[0].isFinal;
            const strongEarthquake = yahoo_data.realTimeData.intensity;
            const regionName = yahoo_data.hypoInfo.items[0].regionName;
            const maxstrongEarthquake = Math.max(...strongEarthquake.split('').map(char => seismicIntensityConversion(char)));
            const calcIntensity = yahoo_data.hypoInfo.items[0].calcintensity.replace("0", "");
            const magnitude = yahoo_data.hypoInfo.items[0].magnitude;
            const depth = yahoo_data.hypoInfo.items[0].depth;
            const Wave_latitude = yahoo_data.psWave.items[0].latitude.replace("N", "");
            const Wave_longitude = yahoo_data.psWave.items[0].longitude.replace("E", "");
            const pRadius = parseFloat(yahoo_data.psWave.items[0].pRadius).toFixed(0);
            const sRadius = parseFloat(yahoo_data.psWave.items[0].sRadius).toFixed(0);
            const report = isFinal ? "final Report" : `Report ${reportNum}`;
            return [currentDateTime, maxstrongEarthquake, apiUrl, report, regionName, calcIntensity, magnitude, depth, Wave_latitude, Wave_longitude, pRadius, sRadius];
        }
    } else {
        console.error('リクエストが失敗しました。');
        return null; // エラー時はnullを返すなど、適切なエラー処理を行う
    }
}