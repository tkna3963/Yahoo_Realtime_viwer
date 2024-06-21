function now_time() {
    currentDateTime = new Date();
    return currentDateTime
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

function Kyoushin_moniter_url() {
    now_time_date = now_time()
    set_time = new Date(now_time_date.getTime() - 3 * 1000);
    url_time_date = Yahoo_Time_date_fromat(set_time)
    const Kyoushin_B_Url = `https://smi.lmoniexp.bosai.go.jp/data/map_img/RealTimeImg/jma_b/${url_time_date}.jma_b.gif`;
    const Tyoushuki_Url = `https://www.lmoni.bosai.go.jp/monitor/data/data/map_img/RealTimeImg/abrspmx_s/${url_time_date}.abrspmx_s.gif`;
    const PGA_Url = `https://smi.lmoniexp.bosai.go.jp/data/map_img/RealTimeImg/acmap_s/${url_time_date}.acmap_s.gif`
    return {
        "Kyoushin_B_Url": Kyoushin_B_Url,
        "Tyoushuki_Url":Tyoushuki_Url,
        "PGA_Url":PGA_Url
    };
}

function Update() {
    const urls = Kyoushin_moniter_url();
    document.getElementById('imgKyoushinB').src = urls.Kyoushin_B_Url;
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

Update();
setInterval(Update, 1000);
updateImages();
