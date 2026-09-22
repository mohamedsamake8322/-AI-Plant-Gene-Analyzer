#!/usr/bin/env python3
"""Small HTTP utilities: session with retries and simple backoff for API collectors.

Usage:
    import request_utils as rq
    resp = rq.get(url, params=params)
"""
from __future__ import annotations

import time
from typing import Optional
from urllib.parse import urlparse

import certifi
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def _build_session(retries: int = 3, backoff_factor: float = 0.5, status_forcelist=(429, 500, 502, 503, 504)) -> requests.Session:
    s = requests.Session()
    s.verify = certifi.where()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=frozenset(["GET", "POST"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    return s


_DEFAULT_SESSION: Optional[requests.Session] = None

# Cadence minimale (en secondes entre deux requêtes) par domaine. KEGG
# bloque explicitement au-delà de 3 requêtes/seconde -- on vise 2/s pour
# garder une marge de sécurité. Les autres domaines gardent le rythme
# par défaut existant.
_MIN_INTERVAL_BY_HOST = {
    "rest.kegg.jp": 0.5,
}
_DEFAULT_MIN_INTERVAL = 0.08
_last_request_time: dict[str, float] = {}


def get_session() -> requests.Session:
    global _DEFAULT_SESSION
    if _DEFAULT_SESSION is None:
        _DEFAULT_SESSION = _build_session()
    return _DEFAULT_SESSION


def _throttle(url: str) -> None:
    host = urlparse(url).netloc
    min_interval = _MIN_INTERVAL_BY_HOST.get(host, _DEFAULT_MIN_INTERVAL)
    last = _last_request_time.get(host, 0.0)
    wait = min_interval - (time.monotonic() - last)
    if wait > 0:
        time.sleep(wait)
    _last_request_time[host] = time.monotonic()


def get(url: str, params: dict | None = None, timeout: int = 30, headers: dict | None = None, retries: int = 3) -> requests.Response:
    """Perform a GET with retries/backoff. Raises requests.RequestException on hard failures."""
    _throttle(url)
    sess = get_session()
    try:
        resp = sess.get(url, params=params, timeout=timeout, headers=headers)
    except Exception:
        delay = 1.0
        for attempt in range(1, retries + 1):
            try:
                resp = requests.get(
                    url, params=params, timeout=timeout, headers=headers,
                    verify=certifi.where(),
                )
                break
            except Exception:
                if attempt == retries:
                    raise
                time.sleep(delay)
                delay *= 2
    resp.raise_for_status()
    return resp