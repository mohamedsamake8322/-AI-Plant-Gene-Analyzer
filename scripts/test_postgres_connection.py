#!/usr/bin/env python3
"""
Diagnostic script to test PostgreSQL connectivity.
"""

from __future__ import annotations

import os
import socket
import sys
from pathlib import Path
from urllib.parse import quote_plus

from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
load_dotenv(ROOT / ".env")

DB_HOST = os.getenv("DB_HOST", "")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "postgres")
DB_USER = os.getenv("DB_USER", "postgres")
DB_PASSWORD = os.getenv("DB_PASSWORD", "")

print("=== PostgreSQL Connectivity Diagnostic ===\n")

print(f"DB_HOST: {DB_HOST}")
print(f"DB_PORT: {DB_PORT}")
print(f"DB_NAME: {DB_NAME}")
print(f"DB_USER: {DB_USER}")
print(f"DB_PASSWORD: {'***' if DB_PASSWORD else '(not set)'}\n")

if not DB_HOST:
    print("❌ ERROR: DB_HOST not configured in .env")
    sys.exit(1)

# Test 1: DNS resolution
print(f"[1] Testing DNS resolution for {DB_HOST}...")
try:
    addr_info = socket.getaddrinfo(DB_HOST, int(DB_PORT))
    print(f"✅ DNS resolved: {addr_info[0][4]}")
except socket.gaierror as e:
    print(f"❌ DNS resolution failed: {e}")
    print("   → Check your internet connection or Supabase domain name")
    sys.exit(1)

# Test 2: Socket connectivity
print(f"\n[2] Testing TCP connection to {DB_HOST}:{DB_PORT}...")
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.settimeout(5)
try:
    sock.connect((DB_HOST, int(DB_PORT)))
    print(f"✅ TCP connection successful")
    sock.close()
except socket.timeout:
    print(f"❌ Connection timeout after 5s")
    print("   → Firewall may be blocking the connection or host is unreachable")
    sys.exit(1)
except socket.error as e:
    print(f"❌ TCP connection failed: {e}")
    print("   → Check firewall rules or if Supabase instance is active")
    sys.exit(1)
finally:
    sock.close()

# Test 3: psycopg connection
print(f"\n[3] Testing psycopg PostgreSQL connection...")
try:
    import psycopg

    database_url = (
        f"postgresql://{quote_plus(DB_USER)}:{quote_plus(DB_PASSWORD)}@{DB_HOST}"
        f":{DB_PORT}/{quote_plus(DB_NAME)}"
    )
    print(f"Connection string: postgresql://{DB_USER}:***@{DB_HOST}:{DB_PORT}/{DB_NAME}")
    
    conn = psycopg.connect(database_url, autocommit=True, connect_timeout=5)
    print(f"✅ PostgreSQL connection successful")
    print(f"   Database version: {conn.info.server_version}")
    conn.close()
except Exception as e:
    print(f"❌ PostgreSQL connection failed: {e}")
    print("   → Check your credentials (user/password) or database name")
    sys.exit(1)

print(f"\n✅ All tests passed! PostgreSQL is accessible.")
