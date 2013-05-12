#!/bin/bash
# Requires bash version 4 or later

getOptions() {
  while getopts "e:t:d:" arg; do
    case $arg in
	t) testname=$OPTARG;;
      	e) exename=$OPTARG;;
      	d) exedir=$OPTARG;;
	?) echo "Error: Too many arguments."
		exit 3;;
    esac
  done
}

evalAndPrint() {
	echo "[/bin/bash]$ $@"
	eval "$@"
}

cd $(dirname $0)

exedir=pwd
getOptions "$@"
testdir=$(pwd)/tests/$testname
inputdir=$testdir/input
outputdir=$testdir/output
accepteddir=$testdir/accepted
result=0

# Ensure a testname was given
if test -z $testname
then
	echo "Error: Name of test not given.  Specify with -t"
	exit 3
fi
# By default, exename is assumed to be same as test name
if test -z $exename
then
	exename=$testname
fi
# Check that exename exists
if test ! -x $exedir/$exename
then
	echo "Error: $exedir/$exename does not exist or has wrong permissions."
	exit 3
fi
# Check that config file exists
if test ! -f $testdir/$testname.conf
then
	echo "Error: Missing config file"
	echo "Must supply $testdir/$testname.conf"
	exit 3
fi


# Source the configuration file
source $testdir/$testname.conf

# Make output directory if it doesn't already exist
if test ! -d $outputdir
then
	mkdir $outputdir
fi

cd $outputdir
if test $? -ne 0 
then
	echo "Unable to create an output directory."
	exit 3
fi

# Clear the output directory
rm -rf *

# Now go to the input directory.  The input directory is the working directory while the test program
# is being executed.
cd $inputdir
if test $? -ne 0 
then
	echo "Input directory does not exist."
	exit 3
fi

# Update input files and accepted files
if test $SVN_UPDATE_FILES = "true" || $SVN_UPDATE_FILES = "TRUE"
then
	svnoptions="--non-interactive --trust-server-cert"
	if test -n "${SVN_REVISION_NUMBER}"
	then
		evalAndPrint svn update $accepteddir $inputdir -r ${SVN_REVISION_NUMBER} $svnoptions
	else
		evalAndPrint svn update $accepteddir $inputdir -r HEAD $svnoptions
	fi
fi

beforefiles=`ls`

# Run program to generate output
evalAndPrint $exedir/$exename ${EXEARGS[@]} 1> /dev/null
exitcode=$?

afterfiles=`ls`

# Move generated files to the output directory
for file1 in ${afterfiles}
do
	matched=1
	for file2 in ${beforefiles}
	do
		if test $file1 = $file2
		then
			matched=0
		fi
	done
	if test ! $matched -eq 0
	then
		evalAndPrint mv $file1 $outputdir
	fi
done

# If program ran successfully, diff output data with accepted data
if test $exitcode -eq 0
then
	echo "${exename^} ran successfully, comparing output with accepted results."

	cd $outputdir

	for(( i=0 ; i<${#ACCEPTED_FILES[@]} ; i++ ))
	do
		echo "*****"
		echo "Checking file ${ACCEPTED_FILES[i]} with ${CMP_PROG[i]}"
		echo "*****"
		# Test that the files are actually present
		if test ! -f $accepteddir/${ACCEPTED_FILES[i]}
		then
			echo "Accepted file ${ACCEPTED_FILES[i]} not present in $accepteddir."
			exit 3
		fi
		if test ! -f ${ACCEPTED_FILES[i]}
		then
			cp $inputdir/${ACCEPTED_FILES[i]} $outputdir
			if test $? -ne 0
			then
				echo "Needed output file ${ACCEPTED_FILES[i]} is not in $outputdir."
				exit 1
			fi
		fi

		# Check that program used for comparison is installed
		which ${CMP_PROG[i]} &> /dev/null
		if test $? -ne 0
		then
			echo "Comparison program, ${CMP_PROG[i]}, is not installed."
			exit 3
		fi

		# Check for differences between the output and accepted file
		evalAndPrint ${CMP_PROG[i]} ${CMP_PROG_ARGS[i]} ${ACCEPTED_FILES[i]} $accepteddir/${ACCEPTED_FILES[i]}
		if test $? -ne 0
		then
			result=1
		fi
	done
else
	echo "${exename^} did not run successfully."
	exit 2
fi

if test $result -eq 0
then
	echo "All output from ${exename^} is consistent with accepted result."
else
	echo "Output from ${exename^} is NOT consistent with accepted result."
fi

exit $result
