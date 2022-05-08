#!/usr/bin/env bash
# Hao, Chunlin@20220426

PIPES=$(dirname $(which metaspectrast))/bin

help_msg () {
	echo "_________ metaSpectraST (v0.0) by Hao, Chunlin _________"
	echo ""
	echo "Usage: metaspectrast cluster [option] -i fragmentation_type <mzML_file>"
	echo "Options:"
	echo "        -i STR    Fragmentation type of the spectra (ETD, HCD, CID-QTOF)"
}

# Check the num of arguments. At least specify one mzML file.
if [[ $# -lt 1 ]]; then help_msg; exit 1; fi

# Loop through arguments
while true; do
	case "$1" in
		-i) fragment=$2; shift 2;;
		-h | --help) help_msg; exit 0; shift 1;;
		*) break;;
	esac
done

# Check fragmentation type
if [[ $fragment != 'HCD' ]] && \
	[[ $fragment != 'ETD' ]] && \
	[[ $fragment != 'CID-QTOF' ]] &&\
	[[ ! -z $fragment ]]; then
	echo "Invalid fragmentation type!"
	echo "Fragmentation type: HCD, ETD, or CID-QTOF"; exit 1;
fi

files=$@
if [[ $# -lt 1 ]]; then 
	echo "Invalid arguments!"; 
	help_msg; exit 1; 
fi

# Find if entire TPP suite is installed
if [ -z "$(which xinteract)" ]; then
	# Import mzML files into SpectraST
	if [[ -z $fragment ]]; then
		for sample in $files; do
			data=$(basename $sample)
			data=${data%.*}
			${PIPES}/spectrast -cn${data} -c_RRS ${sample}
		done
	else
		for sample in $files; do
			data=$(basename $sample)
			data=${data%.*}
			${PIPES}/spectrast -cI${fragment} -cn${data} -c_RRS ${sample}
		done
	fi
	# Combine imported files
	${PIPES}/spectrast -cNimported *.splib
	# Spectral clustering
	${PIPES}/spectrast -cNgrandConsensus -cJU -cAS -c_RRS -c_MGF imported.splib
else
	tpp=$(dirname $(which xinteract))/spectrast
	echo "Invoke [${tpp}]"
	echo ""
	# Import mzML files into SpectraST
	if [[ -z $fragment ]]; then
		for sample in $files; do
			data=$(basename $sample)
			data=${data%.*}
			${tpp} -cn${data} -c_RRS ${sample}
		done
	else
		for sample in $files; do
			data=$(basename $sample)
			data=${data%.*}
			${tpp} -cI${fragment} -cn${data} -c_RRS ${sample}
		done
	fi
	# Combine imported files
	${tpp} -cNimported *.splib
	# Spectral clustering
	${tpp} -cNgrandConsensus -cJU -cAS -c_RRS -c_MGF imported.splib
fi

echo "$(date), job completed!"
