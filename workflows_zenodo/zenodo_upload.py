import os
import sys
import click
import requests
import time

### add deposition info at lines 89-90 (this script works if you already created a deposition draft)
### submit the script with 'python upload_to_zenodo.py /path/to/file_or_folder_to_submit' 
### you will be asked to insert your personal token  

ZENODO_API_URL = "https://zenodo.org/api/deposit/depositions"


class ProgressFileReader:
    """File wrapper to show upload progress without chunking."""

    def __init__(self, file_path):
        self.file = open(file_path, "rb")
        self.file_size = os.path.getsize(file_path)
        self.bytes_read = 0
        self.last_update = time.time()

    def read(self, chunk_size=1024 * 1024):  # read 1 MB at a time internally
        data = self.file.read(chunk_size)
        if not data:
            return b""
        self.bytes_read += len(data)
        now = time.time()
        # Limit print frequency to ~5 updates/sec
        if now - self.last_update > 0.2:
            percent = self.bytes_read / self.file_size * 100
            print(
                f"\rProgress: {percent:.2f}% ({self.bytes_read/1e6:.1f} / {self.file_size/1e6:.1f} MB)",
                end="",
                flush=True,
            )
            self.last_update = now
        return data

    def __len__(self):
        return self.file_size

    def close(self):
        self.file.close()


def upload_file(file_path, deposition_bucket, access_token):
    """
    Upload a file to Zenodo (new API, no manual chunking, with progress bar).
    """
    fname = os.path.basename(file_path)
    file_size = os.path.getsize(file_path)
    url = f"{deposition_bucket}/{fname}"
    headers = {"Authorization": f"Bearer {access_token}"}

    print(f"Uploading {fname} ({file_size / 1e9:.2f} GB)...")

    reader = ProgressFileReader(file_path)
    try:
        response = requests.put(url, headers=headers, data=reader)
    finally:
        reader.close()

    print("\rProgress: 100.00% ✅ Upload complete!")

    if response.status_code not in (200, 201):
        raise RuntimeError(f"❌ Upload failed: {response.status_code}, {response.text}")

    return response.json()


@click.command()
@click.argument(
    "path",
    required=True,
    type=click.Path(exists=True, file_okay=True, dir_okay=True),
)
@click.option("--access-token", prompt=True, help="Zenodo access token")
def main(path, access_token):
    """Upload a file or folder to Zenodo (no chunking, with progress bar)."""
    if path == "-":
        path = sys.stdin.read().strip()

    if not os.path.exists(path):
        click.echo("❌ Provided path does not exist.")
        return

    # Your deposition info
    deposition_id = "<ADD_DEPOSITION_ID>"
    deposition_bucket = "https://zenodo.org/api/files/<ADD_BUCKET_URL>"
    
    critical_size = 50 * 1024 * 1024 * 1024  # 50 GB
    
    if os.path.isfile(path):
        file_size = os.path.getsize(path)
        if file_size > critical_size:
            click.echo(
            f"❌ File size exceeds 50 GB limit (actual: {file_size / 1e9:.2f} GB)"
            )
            return
        try:
            upload_file(path, deposition_bucket, access_token)
        except Exception as e:
            click.echo(str(e))

    elif os.path.isdir(path):
        def get_folder_size(folder_path):
            total_size = 0
            for dirpath, dirnames, filenames in os.walk(folder_path):
                for f in filenames:
                    fp = os.path.join(dirpath, f)
                    total_size += os.path.getsize(fp)
            return total_size

        folder_size = get_folder_size(path)
        if folder_size > critical_size:
            click.echo(
                f"❌ Folder size exceeds 50 GB limit (actual: {folder_size / 1e9:.2f} GB)"
            )
            return
        
        for root, dirs, files in os.walk(path):
            if ".git" in dirs:
                dirs.remove(".git")
            for file in files:
                file_path = os.path.join(root, file)
                try:
                    upload_file(file_path, deposition_bucket, access_token)
                except Exception as e:
                    click.echo(f"\n❌ Failed to upload {file_path}: {e}")
                else:
                    click.echo(f"\n✅ Uploaded {file_path}")
    else:
        click.echo("❌ Path is neither a file nor a folder.")


if __name__ == "__main__":
    main()
