%%20210704
%聂杰文

%%
tic
clear
%%
%-----------------------------接收数据加载---------------------------------
local='D:\OFC2022\OSC_data\OSCdata1012\25m\';
fre=50;
seq=csvread([local 'clear_200M_1ms.csv'],2,3);

Array=seq(:,2);
Tinterval=mean(diff(seq(:,1)));
Length=length(Array);
%%
%--------------------------------信号源设置参数(手动)----------------------------
baud_per_pattern=250;
Fsend=fre/250*1e6;%frequency send :200k
order=2;%PAM 调制阶数
Ttotal=Tinterval*Length;
Npat=Fsend*Ttotal;%number of pattern
sample_per_baud=Length/(baud_per_pattern*Npat);%每个波特内的采样点数
%%
%后期将原始序列分组操作，以求减少积累误差
for i = 1:length(Array)/floor(sample_per_baud)/250
    Array_passage(:,i) = Array(sample_per_baud*250*(i-1)+1:sample_per_baud*250*i);    
end

    %----------------------------------信号采样处理----------------------------
    %等间隔采样，然后硬判决，根据判决结果分组计算方差，以判决线
    %1\二分法将轨道分组，并算出硬判决界限；2、等间隔采样，然后使用之前的硬判决门限将轨道分组,同样是二分法，计算方差，方差之和最小的就是最佳采样点。
    %函数1，输入波形和轨道数，二分法求出各级判决门限；
    %函数2，输入取样之后的波形和各级判决门限，按条件，依次求轨道及方差；求出方差均值

    %----------根据采样后各级门限两旁的方差，寻找最佳采样点（最大方差）--------
for q = 1:length(Array)/floor(sample_per_baud)/250
    Array1 = Array_passage(:,q);
    threshold= Detection(Array1,order);%返回2^n-1个门限值
    for i=1:floor(sample_per_baud)
%         a=Array(round([i:floor(sample_per_baud):end]));
        a=Array1(([i:floor(sample_per_baud):end]));
        variance(i)=variance_calcualtion(a,threshold);
    end
    ind=find(variance==max(variance));
%     Bstream=Array(round([ind:sample_per_baud:end]));%按最佳采样点采样
        Bstream=Array1(([ind:floor(sample_per_baud):end]));%按最佳采样点采样
    if ind>sample_per_baud
        Bstream=[Array1(1); Bstream];
        disp('sample again')
    end
    %----------------比特同步，自相关求峰值------------------------------------
    %%
    source_signal = load('D:\OFC2022\AWG_250bit_PAM4_.csv');
    T_sample=round(length(source_signal)/baud_per_pattern);
    Bstream_source=source_signal(1:T_sample:end);
%     Bstream_source=repmat(Bstream_source,floor(Npat),1);
    %维度矫正
%     Bstream=Bstream(1:length(Bstream_source));

        for shift=1:baud_per_pattern
            Bstream_temp= circshift(Bstream,shift-1);
            R = corrcoef(Bstream_source,Bstream_temp);
            R_result(shift) = R(1,2);
        end
    Shift = find(R_result == max(max(R_result)));
    Bstream_shift= circshift(Bstream,Shift-1);
    diff = xor(Bstream_source,Bstream_shift);
    diff = diff*1;
    sum_q(:,q) = sum(diff);
    %---------------------截尾，转成bit求比特误码率----------------------------
%     Bstream_shift=Bstream_shift(1:baud_per_pattern*(Npat-1));
%     Bit_stream_shift=Baud2bit(Bstream_shift,order,threshold);
%     Bstream_source=Bstream_source(1:baud_per_pattern*(Npat-1));
%     Threshold_source= Detection(Bstream_source,order);%返回2^n-1个门限值
%     Bit_stream_source=Baud2bit(Bstream_source,order,Threshold_source);
% 
%     BER=sum(xor(Bit_stream_source,Bit_stream_shift))/length(Bit_stream_source)
end
BER = sum(sum_q)/Npat/500;