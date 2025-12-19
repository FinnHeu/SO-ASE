# zenodo_upload module

### *class* zenodo_upload.ProgressFileReader(file_path)

Bases: `object`

File wrapper to show upload progress without chunking.

#### close()

#### read(chunk_size=1048576)

### zenodo_upload.upload_file(file_path, deposition_bucket, access_token)

Upload a file to Zenodo (new API, no manual chunking, with progress bar).
