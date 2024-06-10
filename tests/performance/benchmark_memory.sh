for key in 500 1000 1500 2000 2500 3000 3500 4000
do
    # python performance.py benchmark-memory-single --system="unordered" -s $key
    python performance.py benchmark-memory-single --system="ordered" -s $key
done

