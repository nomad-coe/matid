for key in 100 200 300 400 500 600
do
    python performance.py benchmark-memory-single --system="unordered" -s $key
    python performance.py benchmark-memory-single --system="ordered" -s $key
done

