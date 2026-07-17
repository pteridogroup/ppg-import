#!/usr/bin/env python3
"""Download a WorldFerns CSV from Google Drive.

Either a specific file (--file-id) or the newest matching file in a folder
(--folder-id) can be downloaded.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload

WORLDFERNS_NAME_PATTERN = re.compile(r"^WorldFerns_ver_.*\.csv$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--file-id", help="Drive ID of a specific file to download.")
    source.add_argument(
        "--folder-id",
        help=(
            "Drive folder ID to search for the newest file matching "
            "'WorldFerns_ver_*.csv', instead of a specific --file-id."
        ),
    )
    parser.add_argument(
        "--file-name",
        help="Filename to save as; required when using --file-id.",
    )
    parser.add_argument(
        "--out-dir", required=True, help="Directory to save the downloaded file in."
    )
    parser.add_argument("--credentials", required=True)
    args = parser.parse_args()
    if args.file_id and not args.file_name:
        parser.error("--file-name is required when using --file-id")
    return args


def find_latest_file(service, folder_id: str) -> tuple[str, str]:
    query = f"'{folder_id}' in parents and trashed = false"
    res = (
        service.files()
        .list(
            q=query,
            fields="files(id, name, modifiedTime)",
            orderBy="modifiedTime desc",
            pageSize=100,
        )
        .execute()
    )
    matches = [
        f for f in res.get("files", []) if WORLDFERNS_NAME_PATTERN.match(f["name"])
    ]
    if not matches:
        print(
            f"No files matching 'WorldFerns_ver_*.csv' found in folder {folder_id}.",
            file=sys.stderr,
        )
        sys.exit(1)
    latest = matches[0]
    return latest["id"], latest["name"]


def main() -> None:
    args = parse_args()

    scopes = ["https://www.googleapis.com/auth/drive.readonly"]
    creds = service_account.Credentials.from_service_account_file(
        args.credentials,
        scopes=scopes,
    )
    service = build("drive", "v3", credentials=creds)

    if args.folder_id:
        file_id, file_name = find_latest_file(service, args.folder_id)
    else:
        file_id, file_name = args.file_id, args.file_name

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / file_name

    request = service.files().get_media(fileId=file_id)
    with out_path.open("wb") as handle:
        downloader = MediaIoBaseDownload(handle, request)
        done = False
        while not done:
            _, done = downloader.next_chunk()

    print(f"file_id={file_id}")
    print(f"file_name={file_name}")


if __name__ == "__main__":
    main()
