# Name: sra upload
# Author: Emily Lauren Yaklich
# Date: May 07 2025
# Version: Python 3.10
# Description: upload all files in a directory to ncbi sra uploads

import os
from ftplib import FTP
import getpass

# --- CONFIGURATION ---
FTP_HOST = "ftp-private.ncbi.nlm.nih.gov"
FTP_USER = "subftp"
LOCAL_DIR = "/home/emilyyaklich/directory_with_files_to_upload"
REMOTE_DIR = "/uploads/NCBI_username/destination_folder"
EXTENSION = "fastq.gz"
# ----------------------

# Prompt user for password
FTP_PASS = getpass.getpass(prompt="Enter FTP password: ")

# Connect to FTP in active mode
ftp = FTP()
ftp.connect(FTP_HOST, 21)
ftp.login(FTP_USER, FTP_PASS)
ftp.set_pasv(False)

# Change to target directory
ftp.cwd(REMOTE_DIR)

# Upload all matching files
for filename in os.listdir(LOCAL_DIR):
    if filename.endswith(EXTENSION):
        filepath = os.path.join(LOCAL_DIR, filename)
        print(f"Uploading: {filename}")
        with open(filepath, "rb") as f:
            ftp.storbinary(f"STOR {filename}", f)

ftp.quit()
print("All uploads complete.")

