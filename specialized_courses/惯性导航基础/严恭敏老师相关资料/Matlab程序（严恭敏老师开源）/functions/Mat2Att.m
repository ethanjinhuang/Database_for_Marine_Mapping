function att = Mat2Att(M)  
%×ËÌ¬¾ØÕó×ª»¯Îª×ËÌ¬½Ç
	att = [asin(M(3,2)); atan(-M(3,1)/M(3,3)); atan(M(1,2)/M(2,2))] * 180/pi;
    	
    if M(3,3) < 0
        %att(2) = att(2) - sign(att(2))*180;
		if(att(2) < 0)
			att(2) = att(2)+180;
		else
			att(2) = att(2)-180;
        end
    end

	if( M(2,2) > 0)
		if(att(3) < 0)
			att(3) = att(3)+360;
        end
	else
		att(3) = att(3)+180;
    end
