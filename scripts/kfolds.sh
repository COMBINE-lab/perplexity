#!/bin/bash

# set strict mode
set -euo pipefail
IFS=$'\n\t'

OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
left_file=""
right_file=""
output_dir="kfolds"
delete_orig=0
seed=""
k=0
while getopts "h?dl:r:s:o:k:" opt; do 
    case "$opt" in
        h|\?)
            echo "Randomly shuffle and distribute a pair of fasta files into k folds"
            echo "usage: kfolds.sh -l <file> -r <file> -o <directory> -k <int> [-h <file> -f -d -s <string>]"
            echo "IMPORTANT - the output directory will be completely deleted when running!"
            echo "-h    print this help message"
            echo "-l    the left file (1)"
            echo "-r    the right file (2), in the same directory as the left file"
            echo "-k    the number of partitions to make"
            echo "-o    the output directory to save to (default ./kfolds)"
            echo "-f    refresh (overwrite output files)"
            echo "-d    delete input files when done"
            echo "-s    shuffle with a seed (ex: -s 'seed')" # internally through random source and echo
            exit 0
            ;;
        l)  left_file=$OPTARG
            # echo $left_file
            ;;
        r)  right_file=$OPTARG
            # echo $right_file
            ;;
        o)  output_dir=$OPTARG
            ;;
        k)
            k=$OPTARG
            ;;
        d)
            delete_orig=1
            ;;
        s)
            # https://unix.stackexchange.com/questions/496788/does-the-size-of-the-random-source-file-matter
            seed=$OPTARG
            ;;
        esac
done
re='^[0-9]+$'
if ! [[ $k =~ $re ]] || [[ $k -le 1 ]]
then
    echo "k must be a positive integer > 1"
    exit 1
fi

dirname_left=$(dirname left_file)
dirname_right=$(dirname right_file)

if [[ $dirname_left != "$dirname_right" ]]
then
    echo "The left and right reads should be in the same directory"
    exit 1
fi
base_left=$(basename "$left_file")
base_right=$(basename "$right_file")
fl=${base_left%*.fasta}
fr=${base_right%*.fasta}

# to avoid deleting everything in root
#rm -r "${output_dir:?}"
test_dir="$output_dir/test"
train_dir="$output_dir/train"
tmp_dir="$output_dir/tmp"
# a temporary folder

for ((i = 1; i <= k; i++)); do
    mkdir -p "$test_dir/$i"
    mkdir -p "$train_dir/$i"
done
mkdir -p "$tmp_dir"

# https://www.gnu.org/software/coreutils/manual/html_node/Random-sources.html

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

rand_source=''
if [[ -n $seed ]]
then
    # this will be expanded out later
    rand_source='--random-source=<(get_seeded_random $seed)'
fi

# modified from https://superuser.com/questions/396536/how-to-keep-only-every-nth-line-of-a-file
# use an awk script to avoid duplicating the input data twice 

echo "Shuffling and splitting $left_file and $right_file"
carrots=$(grep -c '^>' "$left_file") # keep track of awk progress

cmd="paste $left_file $right_file | paste - - | shuf $rand_source |
    awk -v i=1 -v c=0 -v pd=$tmp_dir -v k=$k -v carrots=$carrots -v percent=0 '{
        print \$1 > pd\"/part\"i\"_1.fasta\";
        print \$3 > pd\"/part\"i\"_1.fasta\";
        print \$2 > pd\"/part\"i\"_2.fasta\";
        print \$4 > pd\"/part\"i\"_2.fasta\";
        i=(i%k)+1;
        c++;
        if (c*100/carrots >= percent) {
            print \"Split \"percent\"% (\"c\"/\"carrots\")\";
            percent += 20;
        }
    }'
"
# all arithmetic in awk is done with floats
# use eval for double string expansion
eval "$cmd"

# copy the partitions into the correct folders
# one extra iteration is not a big deal in the long run
for ((i = 1; i < k+1; i++)); do
    echo "Making fold $i of $k"
    cp "${tmp_dir}/part${i}_1.fasta" "${test_dir}/$i/${fl}_1.fasta"
    cp "${tmp_dir}/part${i}_2.fasta" "${test_dir}/$i/${fr}_2.fasta"
    # find returns a different order every time
    # replace all instances of 1.fasta with 2.fasta to ensure the reads are in the same order
    trainsetone=$(find "$tmp_dir" -type f \( -name '*_1.fasta' -a ! -name "part${i}_1.fasta" \))
    trainsettwo="${trainsetone//1\.fasta/2\.fasta}"
    # confusingly, the > is the redirect *after* 'xargs cat'
    echo "$trainsetone" | xargs cat > "${train_dir}/$i/${fl}_1.fasta"
    echo "$trainsettwo" | xargs cat > "${train_dir}/$i/${fr}_2.fasta"
done

# delete the temporary files
rm -r "${tmp_dir:?}"

if [[ $delete_orig -eq 1 ]]
then
    echo "Deleting original files"
    rm "$left_file" "$right_file"
fi

exit 0