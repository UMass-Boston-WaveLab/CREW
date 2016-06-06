#!/bin/sh
a = RawReceive;
for (( i=1; i<101; i++))
do
echo "Press 1 to run next sample."
echo "Press 2 to exit."
read -n 1 -s -p "Input Selection:" userinput
case $userinput in
    1)
    sleep(60);
    python /home/wavelab/Documents/CREW_EXP_1.py;
    sleep(10)
    n = $(i)
    mv /home/wavelab/Documents/"$a".bin /home/wavelab/Documents/data_1/"$a"_"$n".bin
    ;;
    2)
    echo "Closing script."
    break
    ;;
esac
done
