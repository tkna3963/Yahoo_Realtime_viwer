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
            "0": "Required_files\\Shindo_0.png",
            "1": "Required_files\\Shindo_1.png",
            "2": "Required_files\\Shindo_2.png",
            "3": "Required_files\\Shindo_3.png",
            "4": "Required_files\\Shindo_4.png",
            "5-": "Required_files\\Shindo_5-.png",
            "5+": "Required_files\\Shindo_5+.png",
            "6-": "Required_files\\Shindo_6-.png",
            "6+": "Required_files\\Shindo_6+.png",
            "7": "Required_files\\Shindo_7.png"
        };

        return levelMap[String(scaleData)] || "#0d0d66";
    } catch (error) {
        console.error(error);
    }
}

function seismicIntensityscaleColorConversion(scaleData) {
    try {
        console.log(scaleData)
        const colorMap = {
            null: "#0d0d66",
            "-1": "#9b9b1b",
            "0": "#0fb02b",
            "0.5": "#f4e200",
            "1": "#fbc300",
            "1.5": "#ffaf00",
            "2": "#f90",
            "2.5": "#ff7e00",
            "3": "#ff6200",
            "3.5": "#fc4c02",
            "4": "#f53605",
            "4.5": "#f11520",
            "5": "#ed0047",
            "5.5": "#e30071",
            "6": "#dc009c",
            "6.5": "#c900ba",
            "7": "#b600d7"
        };
        var result_data=colorMap[String(scaleData)] || "#0d0d66";
        return result_data
    } catch (error) {
        console.error(error);
    }
}

function calculateZoomLevel(distance) {
    // 距離に基づいて適切なズーム レベルを計算する
    if (distance <= 100) return 10;   // 距離が 100km 未満の場合、ズーム レベル 10 を返す
    if (distance <= 500) return 8;    // 距離が 500km 未満の場合、ズーム レベル 8 を返す
    if (distance <= 1000) return 6;   // 距離が 1000km 未満の場合、ズーム レベル 6 を返す
    if (distance <= 1500) return 5;   // 距離が 1500km 未満の場合、ズーム レベル 5 を返す
    return 5;                         // それ以外の場合、ズーム レベル 4 を返す
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

function formatTimeDate(dateString, timeZone = 'Asia/Tokyo') {
    const options = { year: 'numeric', month: '2-digit', day: '2-digit', hour: '2-digit', minute: '2-digit', second: '2-digit', timeZone: timeZone };
    const date = new Date(dateString);
    const formatter = new Intl.DateTimeFormat('ja-JP', options);
    const formattedDate = formatter.format(date);
    return formattedDate.replace(/\//g, '年').replace(' ', '日').replace(':', '時').replace(':', '分') + '秒';
}

let set_time_counter = 0;
function yahooShingenn() {
    var timeMachineInput = document.getElementById("time_machine");
    var timeMachineValue = timeMachineInput.value;

    var set_time;
    if (timeMachineValue === "") {
        currentDateTime = new Date();
        set_time = new Date(currentDateTime.getTime() - 3 * 1000);
    } else {
        set_time = new Date(timeMachineValue);
        set_time_counter += 1;
        set_time = new Date(set_time.getTime() + set_time_counter * 1000);
        console.log(set_time, set_time_counter)
    }

    const time_set = formatDateTimeForUrl(set_time);
    const apiUrl = `https://weather-kyoshin.west.edge.storage-yahoo.jp/RealTimeData/${time_set}.json`;
    const xhr = new XMLHttpRequest();
    xhr.open("GET", apiUrl, false);
    xhr.send();

    if (xhr.readyState === 4 && xhr.status === 200) {
        const yahoo_data = JSON.parse(xhr.responseText);
        if (yahoo_data.hypoInfo === null) {
            const strongEarthquake = yahoo_data.realTimeData.intensity;
            const maxstrongEarthquake = Math.max(...strongEarthquake.split('').map(char => seismicIntensityConversion(char)));
            const result_data = [set_time, maxstrongEarthquake, apiUrl, "緊急地震速報は出ていません。", 0, 0, 0, 0, 36.0, 137.9, 0, 0];
            console.log(result_data)
            return result_data
        } else {
            const reportId=yahoo_data.hypoInfo.items[0].reportId;
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
            let pRadius = parseFloat(yahoo_data.psWave.items[0].pRadius).toFixed(0);
            let sRadius = parseFloat(yahoo_data.psWave.items[0].sRadius).toFixed(0);

            if (isNaN(pRadius)) {
                pRadius = 2000;
            }
            if (isNaN(sRadius)) {
                sRadius = 2000;
            }

            let report;
            if (isFinal === "true") {
                report = `${reportId} 第${reportNum}報(最終報)`;
            } else {
                report = `${reportId} 第${reportNum}報`;
            }
            result_data = [set_time, maxstrongEarthquake, apiUrl, report, regionName, calcIntensity, magnitude, depth, Wave_latitude, Wave_longitude, pRadius, sRadius];
            console.log(result_data)
            return result_data
        }
    } else {
        console.error('リクエストが失敗しました。');
        return null; // エラー時はnullを返すなど、適切なエラー処理を行う
    }
}

function Loc_data_reader() {
    // XMLHttpRequestを使用して同期的にJSONファイルを読み込む
    const xhr = new XMLHttpRequest();
    xhr.open('GET', 'https://weather-kyoshin.west.edge.storage-yahoo.jp/SiteList/sitelist.json', false); // 同期的にリクエストを行う
    xhr.send();
    json_data=JSON.parse(xhr.responseText);
    if (xhr.status === 200) {
        return json_data; // JSONを解析して返す
    } else {
        console.error('Error loading settings:', xhr.status);
        return null;
    }
}

function initAndUpateChart(newData) {
    var ctx = document.getElementById('realtime-chart').getContext('2d');
    if (!window.chart) {
        window.chart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: 'リアルタイム強震最大震度',
                    data: [],
                    borderColor: seismicIntensityscaleColorConversion(newData["y"]),
                    borderWidth: 2,
                    pointRadius: 0
                }]
            },
            options: {
                scales: {
                    y: {
                        beginAtZero: true,
                        suggestedMax: 7
                    }
                }
            }
        });
    }

    // 最新データを追加
    window.chart.data.labels.push(newData.x.toLocaleTimeString());
    window.chart.data.datasets[0].data.push(newData.y);
    
    // チャートを更新
    window.chart.update();
}
