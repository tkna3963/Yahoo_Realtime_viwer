User
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

function scaleColorConversion(scaleData) {
    try {
        const colorMap = {
            null: "#0d0d66",
            "-1": "#9b9b1b",
            "0": "#0fb02b",
            "1": "#fbc300",
            "2": "#f90",
            "3": "#ff6200",
            "4": "#f53605",
            "5-": "#ed0047",
            "5+": "#e30071",
            "6-": "#dc009c",
            "6+": "#c900ba",
            "7": "#b600d7"
        };
        return colorMap[String(scaleData)] || "#0d0d66";
    } catch (error) {
        console.error(error);
    }
}

function scaleLevalConversion(scaleData) {
    try {
        const levelMap = {
            null: "#0d0d66",
            "0": "人は揺れを感じないが、地震計には記録される。",
            "1": "屋内で静かにしている人の中には、揺れをわずかに感じる人がいる。",
            "2": "屋内で静かにしている人の大半が、揺れを感じる。眠っている人の中には、目を覚ます人もいる｡",
            "3": "	屋内にいる人のほとんどが、揺れを感じる。歩いている人の中には、揺れを感じる人もいる。眠っている人の大半が、目を覚ます。",
            "4": "ほとんどの人が驚く。歩いている人のほとんどが、揺れを感じる。眠っている人のほとんどが、目を覚ます。",
            "5-": "大半の人が、恐怖を覚え、物につかまりたいと感じる。",
            "5+": "大半の人が、物につかまらないと歩くことが難しいなど、行動に支障を感じる。",
            "6-": "立っていることが困難になる。",
            "6+": "立っていることができず、はわないと動くことができない。揺れにほんろうされ、動くこともできず、飛ばされることもある。",
            "7": "立っていることができず、はわないと動くことができない。揺れにほんろうされ、動くこともできず、飛ばされることもある。"
        };

        return levelMap[String(scaleData)] || "#0d0d66";
    } catch (error) {
        console.error(error);
    }
}

function scaleLevalImgConversion(scaleData) {
    try {
        const levelMap = {
            "0": "Shindo_img\\Shindo_0.png",
            "1": "Shindo_img\\Shindo_1.png",
            "2": "Shindo_img\\Shindo_2.png",
            "3": "Shindo_img\\Shindo_3.png",
            "4": "Shindo_img\\Shindo_4.png",
            "5-": "Shindo_img\\Shindo_5-.png",
            "5+": "Shindo_img\\Shindo_5+.png",
            "6-": "Shindo_img\\Shindo_6-.png",
            "6+": "Shindo_img\\Shindo_6+.png",
            "7": "Shindo_img\\Shindo_7.png"
        };

        return levelMap[String(scaleData)] || "#0d0d66";
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
    const currentDateTime = new Date();
    //const currentDateTime = new Date(2024, 0, 1, 16,11,0);
    const fiveSecondsAgo = new Date(currentDateTime.getTime() - 0 * 1000);
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
            const result_data=[currentDateTime, maxstrongEarthquake, apiUrl, "緊急地震速報は出ていません。", 0, 0, 0, 0, 36.0, 137.9, 0, 0];
            console.log(result_data)
            return result_data
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
            let report;
            if (isFinal === true) {
                report = `最終報(第${reportNum}報)`;
            } else {
                report = `第${reportNum}報`;
            }
            result_data=[currentDateTime, maxstrongEarthquake, apiUrl, report, regionName, calcIntensity, magnitude, depth, Wave_latitude, Wave_longitude, pRadius, sRadius];
            console.log(result_data)
            return result_data
        }
    } else {
        console.error('リクエストが失敗しました。');
        return null; // エラー時はnullを返すなど、適切なエラー処理を行う
    }
}

function toggleElements() {
    var earthquakeInfo = document.getElementById("earthquake_info");
    if (earthquakeInfo.style.display === "none") {
        earthquakeInfo.style.display = "block";
    } else {
        earthquakeInfo.style.display = "none";
    }
}

var isElementsVisible = true;

function toggleElements() {
    var earthquakeInfo = document.getElementById("earthquake_info");
    if (isElementsVisible) {
        earthquakeInfo.style.display = "none";
    } else {
        earthquakeInfo.style.display = "block";
    }
    isElementsVisible = !isElementsVisible;
}