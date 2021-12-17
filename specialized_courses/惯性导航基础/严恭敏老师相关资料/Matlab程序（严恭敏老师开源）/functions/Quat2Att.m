function Att = Quat2Att(q)  
%姿态四元数转化为姿态角
    Att = Mat2Att(Quat2Mat(q));
