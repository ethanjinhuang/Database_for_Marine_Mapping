function qnb = Att2Quat(att)  
%姿态角转化为姿态四元数
    qnb = Mat2Quat(Att2Mat(att));
