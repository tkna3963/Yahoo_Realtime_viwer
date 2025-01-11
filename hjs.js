function updateCircleInfo(data, details) {
    document.getElementById("max-intensity").textContent = data.MAXseismicIntensity;
    if (data.yahoo_datas.situation === "EEW_has_been_issued") {
        document.getElementById("calc-intensity").textContent = data.yahoo_datas.calcIntensity || "--";
        document.getElementById("magnitude").textContent = data.yahoo_datas.magnitude || "--";
        document.getElementById("p-radius").textContent = `${data.yahoo_datas.pRadius}km(約${(Math.round(data.SFCD / 8)) - calculateTimeDifference(parseCustomDate(data["yahoo_datas"]["get_date"]), parseCustomDate(details.originTime))}秒)` || "--";
        document.getElementById("s-radius").textContent = `${data.yahoo_datas.sRadius}km(約${(Math.round(data.SFCD / 4)) - calculateTimeDifference(parseCustomDate(data["yahoo_datas"]["get_date"]), parseCustomDate(details.originTime))}秒)` || "--";
        document.getElementById("reportNum").textContent = data.yahoo_datas.reportNum || "--";
        document.getElementById("calclocate-intensity").textContent = data.SFCI || "--";
    } else {
        document.getElementById("calc-intensity").textContent = "--";
        document.getElementById("magnitude").textContent = "--";
        document.getElementById("p-radius").textContent = "--";
        document.getElementById("s-radius").textContent = "--";
        document.getElementById("reportNum").textContent = "--";
        document.getElementById("calclocate-intensity").textContent = "--";
    }
}

function main() {
    const datas_bords = datas_bord();

    document.getElementById("get_time").innerText = datas_bords["yahoo_datas"]["get_date"];

    const latitude = datas_bords.latitude || "不明";
    const longitude = datas_bords.longitude || "不明";
    document.getElementById("location").innerText = `緯度: ${latitude}, 経度: ${longitude}`;

    const detailsContainer = document.getElementById("earthquake-details");
    detailsContainer.innerHTML = "";

    let details = null;

    if (datas_bords.yahoo_datas.situation === "EEW_has_been_issued") {
        details = datas_bords.yahoo_datas;

        const detailItems = [
            { label: "報告ID", value: details.reportId },
            { label: "発報番号", value: details.reportNum },
            { label: "報告時間", value: details.reportTime },
            { label: "地震発生時刻", value: details.originTime },
            { label: "震源地名", value: details.regionName },
            { label: "緯度経度座標", value: `${details.Wave_latitude},${details.Wave_longitude}` },
            { label: "震源距離", value: `${datas_bords.SFCD}km` },
            { label: "震源の深さ (km)", value: details.depth },
            { label: "予想最大深度", value: details.calcIntensity },
            { label: "マグニチュード", value: details.magnitude },
            { label: "波の地点", value: `${details.pRadius}km(約${(Math.round(datas_bords.SFCD / 8)) - calculateTimeDifference(parseCustomDate(datas_bords["yahoo_datas"]["get_date"]), parseCustomDate(details.originTime))}秒) ${details.sRadius}km(約${(Math.round(datas_bords.SFCD / 4)) - calculateTimeDifference(parseCustomDate(datas_bords["yahoo_datas"]["get_date"]), parseCustomDate(details.originTime))}秒)` }
        ];

        detailItems.forEach(item => {
            const div = document.createElement("div");
            div.className = "data-point";
            div.innerHTML = `<label>${item.label}:</label><value>${item.value}</value>`;
            detailsContainer.appendChild(div);
        });
    } else {
        detailsContainer.innerHTML = "<div class='data-point'>現在有効な地震情報はありません。</div>";
    }

    document.getElementById("jsonarea").value = objectToJsonText(datas_bords);
    updateCircleInfo(datas_bords, details);
}

setInterval(main, 1000);

