rm -f CRF/crf_test; echo "rm CRF/crf_test"
rm -f CRF/crf_test.o; echo "rm CRF/crf_test.o"
cd CRF; echo "cd CRF" 
echo "make crf_test;"
make crf_test; echo "cd .."
cd ..; echo "done"
