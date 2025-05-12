#!/bin/bash
# -- baseFolder
baseFolder="$PWD"
# -- listOfFiles
listOfFiles="$PWD/$1"
# -- tree name
treeName="jetTree"
# -- root macro
rootMacro="runPicoHFJetMaker.C"
# -- bad run list file
badRunListFileName="BadRunList_14.list"
# -- production Id
productionId=$(date +%F_%H-%M)
# -- set STAR software version
starVersion="pro"
# --max number of files
maxNFiles=1

#================================================================
# -- submission xml file
templateXml="template.xml"

# -- job submission directory
mkdir -p "${baseFolder}/submit/${productionId}/job"

# -- result directory
mkdir -p "${baseFolder}/submit/${productionId}/production"

cd "${baseFolder}/submit/${productionId}/job"

# -- prepare folder
mkdir -p report err log list csh

# -- check for prerequisites and create links
folders=(".sl73_gcc485" "StRoot")

echo -n "Checking prerequisites folders ...  "
for folder in "${folders[@]}"; do
    if [ ! -d "${baseFolder}/${folder}" ]; then
        echo "${folder} does not exist in ${baseFolder}"
        exit 1
    else
        ln -sf "${baseFolder}/${folder}"
    fi
done
echo "ok"

# -- check for run macro
echo -n "Checking run macro ...             "
if [ ! -e "${baseFolder}/StRoot/macros/${rootMacro}" ]; then
    echo "${rootMacro} does not exist in ${baseFolder}/StRoot/macros"
    exit 1
fi
echo "ok"

# -- check for xml file
echo -n "Checking xml file  ...             "
if [ ! -e "${baseFolder}/submit/${templateXml}" ]; then
    echo "XML ${templateXml} does not exist"
    exit 1
else
    ln -sf "${baseFolder}/submit/${templateXml}"
fi
echo "ok"

# -- check for bad run list
echo -n "Checking bad run list  ...         "
if [ -e "${baseFolder}/${badRunListFileName}" ]; then
    cp "${baseFolder}/${badRunListFileName}" BadRunList_14.list
elif [ -e "${baseFolder}/picoLists/${badRunListFileName}" ]; then
    cp "${baseFolder}/picoLists/${badRunListFileName}" BadRunList_14.list
else
    echo "${badRunListFileName} does not exist in ${baseFolder} nor ${baseFolder}/picoLists"
    exit 1
fi
echo "ok"

# -- check listOfFiles file list
echo -n "Checking listOfFiles file list ...       "
if [ ! -e "${listOfFiles}" ]; then
    echo "Filelist ${listOfFiles} does not exist"
    exit 1
fi

# -- cleanup old local packages
[ -e LocalLibraries.zip ] && rm LocalLibraries.zip
[ -d LocalLibraries.package ] && rm -rf LocalLibraries.package

# -- submit
generatedXml="generated.xml"
[ -e "${generatedXml}" ] && rm "${generatedXml}"

cat <<EOF >"${generatedXml}"
<?xml version="1.0" encoding="utf-8" ?>
<!DOCTYPE note [
<!ENTITY rootMacro "${rootMacro}">
<!ENTITY treeName "${treeName}">
<!ENTITY productionId "${productionId}">
<!ENTITY baseFolder "${baseFolder}">
<!ENTITY badRunListFileName "${badRunListFileName}">
<!ENTITY listOfFiles "${listOfFiles}">
<!ENTITY starVersion "${starVersion}">
<!ENTITY maxNFiles "${maxNFiles}">
]>
EOF
# -- add the rest of the xml file except the first line <?xml version="1.0" encoding="utf-8" ?>
tail -n +2 "${templateXml}" >>"${generatedXml}"

star-submit-beta "${generatedXml}"
