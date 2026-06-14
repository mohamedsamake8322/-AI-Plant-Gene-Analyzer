#!/usr/bin/env python3
"""Small HTTP utilities: session with retries and simple backoff for API collectors.

Usage:
    import request_utils as rq
    resp = rq.get(url, params=params)
"""
from __future__ import annotations

import time
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


def _build_session(retries: int = 3, backoff_factor: float = 0.5, status_forcelist=(429, 500, 502, 503, 504)) -> requests.Session:
    s = requests.Session()
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


def get_session() -> requests.Session:
    global _DEFAULT_SESSION
    if _DEFAULT_SESSION is None:
        _DEFAULT_SESSION = _build_session()
    return _DEFAULT_SESSION


def get(url: str, params: dict | None = None, timeout: int = 30, headers: dict | None = None, retries: int = 3) -> requests.Response:
    """Perform a GET with retries/backoff. Raises requests.RequestException on hard failures.

    The underlying urllib3 Retry will handle retries for 429/5xx. We also apply a small sleep
    after each successful request to be polite to public APIs.
    """
    sess = get_session()
    try:
        resp = sess.get(url, params=params, timeout=timeout, headers=headers)
    except Exception:
        # Last-resort: fall back to a plain requests.get with simple exponential backoff
        delay = 1.0
        for attempt in range(1, retries + 1):
            try:
                resp = requests.get(url, params=params, timeout=timeout, headers=headers)
                break
            except Exception:
                if attempt == retries:
                    raise
                time.sleep(delay)
                delay *= 2
    # If server replies 429, 5xx etc., let callers decide; session retried already.
    # Small polite pause to reduce burst rate
    try:
        time.sleep(0.08)
    except Exception:
        pass
    resp.raise_for_status()
    return resp
