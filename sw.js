const CACHE_NAME = 'yr-pwa-v1';
const urlsToCache = [
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/index.html',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/script.js',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/styles.css',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/manifest.json',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/Required_files/epicenter.png',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/Required_files/centerARVs.json',
    'https://tkna3963.github.io/Yahoo_Realtime_viwer/Required_files/japanmap.json'
];

self.addEventListener('install', event => {
    event.waitUntil(
        caches.open(CACHE_NAME)
            .then(cache => cache.addAll(urlsToCache))
    );
});

self.addEventListener('fetch', event => {
    event.respondWith(
        caches.match(event.request)
            .then(resp => resp || fetch(event.request))
    );
});
