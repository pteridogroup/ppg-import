#!/usr/bin/env python3
"""Download a WorldFerns CSV from Google Drive by file ID."""

from __future__ import annotations

import argparse
from pathlib import Path

from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--file-id", required=True)
    parser.add_argument("--out-path", required=True)
    parser.add_argument("--credentials", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    scopes = ["https://www.googleapis.com/auth/drive.readonly"]
    creds = service_account.Credentials.from_service_account_file(
        args.credentials,
        scopes=scopes,
    )
    service = build("drive", "v3", credentials=creds)

    out_path = Path(args.out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    request = service.files().get_media(fileId=args.file_id)
    with out_path.open("wb") as handle:
        downloader = MediaIoBaseDownload(handle, request)
        done = False
        while not done:
            _, done = downloader.next_chunk()

    print(out_path.name)


if __name__ == "__main__":
    main()
