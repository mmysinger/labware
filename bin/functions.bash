
f_head()
{
	echo ${1%/*}
}

f_tail()
{
	echo ${1##*/}
}

f_root()
{
	echo ${1%.*}
}

f_ext()
{
	echo ${1##*.}
}

