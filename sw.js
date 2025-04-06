let hostOrigin = '';

const CACHE_NAME = 'yr-pwa-v1';
let urlsToCache = [
    '/', 
    '/index.html',
    '/scripts.js',
    '/styles.css',
    '/manifest.json',
    '/Required_files/epicenter.png',
    '/Required_files/centerARVs.json',
    '/Required_files/japanmap.json'
];

// ホスト情報を受け取る
self.addEventListener('message', event => {
    if (event.data.type === 'SET_HOST') {
        hostOrigin = event.data.host;
        // 必要なら hostOrigin を使って URL を補完したりログに出力したりできる
        console.log('Received host:', hostOrigin);
    }
});

self.addEventListener('install', event => {
    event.waitUntil(
        caches.open(CACHE_NAME).then(cache => {
            return cache.addAll(urlsToCache.map(url => {
                // フルパスにしたい場合はこちら
                return (url.startsWith('http') ? url : hostOrigin + url);
            }));
        })
    );
});

self.addEventListener('fetch', event => {
    event.respondWith(
        caches.match(event.request).then(response => {
            return response || fetch(event.request);
        })
    );
});
