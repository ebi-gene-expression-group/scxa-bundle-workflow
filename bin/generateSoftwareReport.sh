#!/usr/bin/env bash

subworkflow=$1
outfile=$2

template=${NXF_ASSETS}/${NXF_ORG}/$subworkflow/conf/software.tsv

if [ ! -e "$template" ]; then
    echo "Software versions file $template does not exist" 1>&2
    exit 1
fi

# The software config file lists the software, environments and packages used.
# We can then extract the versions from the environments and create a properly
# versioned file for the SCXA interface

echo -e "Analysis\tSoftware\tVersion\tCitation" > $outfile
tail -n +2 $template | while read -r l; do
    analysis=$(echo "$l" | awk -F'\t' '{print $1}')
    software=$(echo "$l" | awk -F'\t' '{print $2}')
    environment=$(echo "$l" | awk -F'\t' '{print $3}')
    package=$(echo "$l" | awk -F'\t' '{print $4}')
    citation=$(echo "$l" | awk -F'\t' '{print $5}')

    environment_file=${NXF_ASSETS}/${NXF_ORG}/$subworkflow/envs/${environment}.yml
    if [ ! -e "$environment_file" ]; then
        echo "Environment file $environment_file not found" 1>&2
        exit 1
    fi

    version_line=$(grep -ri "${package}=" $environment_file)
    if [ $? -ne 0 ]; then
        echo "No version found for package $package in environment file $environment_file" 1>&2
        exit 1
    fi
    version=$(echo $version_line | sed 's/\- //' | awk -F'=' '{print $2}')

    echo -e "$analysis\t$software\t$version\t$citation" >> $outfile   
done

