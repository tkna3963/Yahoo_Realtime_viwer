const CACHE_NAME = 'yr-pwa-v1';
const urlsToCache = [
    '/',
    '/index.html',
    '/scripts.js',
    '/styles.css',
    '/manifest.json',
    '/Required_files/epicenter.png',
    '/Required_files/centerARVs.json',
    '/Required_files/japanmap.json'
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
