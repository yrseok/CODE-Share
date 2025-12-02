#!/bin/bash
#===========< Before run, you need to confirm API on Website >
# Written by Yerim Seok #20240731

#======================< Edit Only Here >==========================
#authKey	use your own key
#key="S8Lzd4M3QUKC83eDN-FCVQ"	#yrs

# 시작 날짜와 시간 설정
start_date="202104181600"
end_date="202104181600"

# 수집 주기 설정 (current setting: 1시간 간격-for hourly data)
#calculate date well
#min data -> hours=0 & adjust minute    2, 10, 30, ...
step_day=0
step_hours=1
step_minute=0

#save directory path
save_path="/home2/GK2A/L1B/IR133/FD/"
#file type LE1B LE2 VI004 EA KO etc
#change file name "GK2A/LE1B/WV063/FD/"
file_name="GK2A/LE1B/IR133/FD/"
save_file="gk2a_ami_le1b_ir133_fd020ge_"
#========================== < Do not Edit > =============================
# 현재 날짜와 시간 초기화
current_time=$start_date

while [ "$current_time" -le "$end_date" ]; do

	OUT=$save_path$save_file${current_time}".nc"    #file check
	echo $OUT


  if [ -e "$OUT" ]; then	# file O

        echo "$current_time" "File exists skipping download."

 else		# file X
  
  echo "$current_time" "File does not exist. Now start download."
	
  wget --no-check-certificate --user-agent Mozilla/4.0 --content-disposition --ignore-length -P "$save_path" "https://apihub-pub.kma.go.kr/api/typ05/api/${file_name}data?date=${current_time}&authKey=$key" || true


 fi		#end if

 # 다음 시간 계산
  current_time=`date -d "${current_time:0:4}-${current_time:4:2}-${current_time:6:2} ${current_time:8:2}:${current_time:10:2}:00 ${step_day} day ${step_hours} hours ${step_minute} minute" +"%Y%m%d%H%M"`

done
