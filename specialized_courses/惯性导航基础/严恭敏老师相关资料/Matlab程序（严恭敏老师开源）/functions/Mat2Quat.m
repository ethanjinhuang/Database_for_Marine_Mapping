function qnb = Mat2Quat(M)  
%姿态矩阵转化为姿态四元数   
    qnb = [                     sqrt(abs(1.0 + M(1,1) + M(2,2) + M(3,3)))/2.0;
          sign(M(3,2)-M(2,3)) * sqrt(abs(1.0 + M(1,1) - M(2,2) - M(3,3)))/2.0;
          sign(M(1,3)-M(3,1)) * sqrt(abs(1.0 - M(1,1) + M(2,2) - M(3,3)))/2.0;
          sign(M(2,1)-M(1,2)) * sqrt(abs(1.0 - M(1,1) - M(2,2) + M(3,3)))/2.0 ]; 
