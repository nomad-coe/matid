for key in 500 1000 1500 2000 2500 3000 3500 4000
do
    python performance.py benchmark-cpu-single --system="ordered" -s $key
done
for key in 100 200 300 400 500
do
    python performance.py benchmark-cpu-single --system="unordered" -s $key
done

