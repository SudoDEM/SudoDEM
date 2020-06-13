for file in ./*; do
    fname="$(basename "$file")"
    prefix="${fname:0:3}"
    if [ "lib" = $prefix ]; then
	if [[ $fname = *".so"* ]]; then #it is a dynamic library
            echo $fname
	    patchelf --set-rpath '$ORIGIN' $fname
	fi
    fi
done

cd ./py/PyQt4
for file in ./*; do
    fname="$(basename "$file")"
    prefix="${fname:0:2}"
    if [ "Qt" = $prefix ]; then
	if [[ $fname = *".so"* ]]; then #it is a dynamic library
            echo $fname
	    patchelf --set-rpath '$ORIGIN/../../' $fname
	fi
    fi
done
