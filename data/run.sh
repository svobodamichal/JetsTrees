#!/bin/bash
productionId=$(date +%F)
# delete old production directory if exists

if [[ -d "submit/${productionId}" ]]; then
    echo "deleting old production directory"
    rm -rf "submit/${productionId}"
fi

filelists=(
    "filelists/embedding/pt3_5.list"
    "filelists/embedding/pt5_7.list"
    "filelists/embedding/pt7_9.list"
    "filelists/embedding/pt9_11.list"
    "filelists/embedding/pt11_15.list"
    "filelists/embedding/pt15_20.list"
    "filelists/embedding/pt20_25.list"
    "filelists/embedding/pt25_30.list"
    "filelists/embedding/pt30_40.list"
    "filelists/embedding/pt40_50.list"
    "filelists/embedding/pt50_-1.list"
)
for list_file in "${filelists[@]}"; do
    ./submit/submit.sh "$list_file"
    if [[ $? -ne 0 ]]; then
        echo "Error processing $list_file"
        exit 1
    fi
done
