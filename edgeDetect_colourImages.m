clear all;
%file1 = input('Enter input file name:','s');
file1 = 'pic(3).jpg';
obj = imread(file1);
image = obj;
[r, c, d] = size(obj);
alph = 1;
timeTaken = cputime;
thresh = 30;

%normalizing factors
s_sym = ((1-exp(-alph))^2)/(1+2*alph*exp(-alph)-exp(-2*alph));
s_anti = (1-exp(-alph))*(1-exp(-alph))/(exp(-alph));

for i=1:r
    
    yNeg = zeros(c,d);
    yPos = zeros(c,d);
    for d1 = 1:d
        yNeg(1,d1) = obj(i,1,d1)/2 + s_sym*obj(i,1,d1)/2;
        yNeg(2,d1) = obj(i,2,d1)/2 + s_sym*obj(i,2,d1)/2;
        yPos(c,d1) = obj(i,c,d1)/2 + s_sym*obj(i,c,d1)/2;
        yPos(c-1,d1) = obj(i,c-1,d1)/2 + s_sym*obj(i,c-1,d1)/2;
    end
    
    for n = 3:c
        for d1 = 1:d
            yNeg(n,d1) = s_sym*obj(i,n,d1) - exp(-alph)*(1-alph)*s_sym*obj(i,n-1,d1) + 2*exp(-alph)*yNeg(n-1,d1)-exp(-2*alph)*yNeg(n-2,d1);
        end
    end
    
    for n = c-2:-1:1
        for d1 = 1:d
            yPos(n,d1) = s_sym*obj(i,n,d1) - exp(-alph)*(1-alph)*s_sym*obj(i,n+1,d1) + 2*exp(-alph)*yPos(n+1,d1)-exp(-2*alph)*yPos(n+2,d1);
        end
    end
    for n = 1:c
        for d1 = 1:d
            obj1(i,n,d1) = yNeg(n,d1) + yPos(n,d1) - s_sym*obj(i,n,d1);
        end
    end
    
    clear yNeg;
    clear yPos;
end

obj = obj1;
clear obj1;

for j = 1:c
    
    yNeg = zeros(r,d);
    yPos = zeros(r,d);
    
    for n = 3:r
        for d1 = 1:d
            yNeg(n, d1) = s_anti * exp(-alph) * obj(n-1,j,d1) + 2 * exp(-alph) * yNeg(n-1, d1) - exp(-2 * alph) * yNeg(n-2, d1) ;
        end
    end
    for n = r-2:-1:1
        for d1 = 1:d
            yPos(n, d1) = s_anti * exp(-alph) * obj(n+1,j,d1) + 2 * exp(-alph) * yPos(n+1, d1) - exp(-2 * alph) * yPos(n+2, d1) ;
        end
    end
    
    for n = 1:r
        for d1 = 1:d
            obj1(n,j,d1) = yNeg(n,d1) - yPos(n,d1);
        end
    end
    
end

Iy = obj1;
obj = image;
clear obj1;

for j=1:c
    
    yNeg = zeros(r,d);
    yPos = zeros(r,d);
    for d1 = 1:d
        yNeg(1,d1) = obj(i,1,d1)/2 + s_sym*obj(i,1,d1)/2;
        yNeg(2,d1) = obj(i,2,d1)/2 + s_sym*obj(i,2,d1)/2;
        yPos(c,d1) = obj(i,c,d1)/2 + s_sym*obj(i,c,d1)/2;
        yPos(c-1,d1) = obj(i,c-1,d1)/2 + s_sym*obj(i,c-1,d1)/2;
    end
    
    for n = 3:r
        for d1 = 1:d
            yNeg(n,d1) = s_sym*obj(n,j,d1) - exp(-alph)*(1-alph)*s_sym*obj(n-1,j,d1) + 2*exp(-alph)*yNeg(n-1,d1)-exp(-2*alph)*yNeg(n-2,d1);
        end
    end
    
    for n = r-2:-1:1
        for d1 = 1:d
            yPos(n,d1) = s_sym*obj(n,j,d1) - exp(-alph)*(1-alph)*s_sym*obj(n+1,j,d1) + 2*exp(-alph)*yPos(n+1,d1)-exp(-2*alph)*yPos(n+2,d1);
        end
    end
    for n = 1:r
        for d1 = 1:d
            obj1(n,j,d1) = yNeg(n,d1) + yPos(n,d1) - s_sym*obj(n,j,d1);
        end
    end
    
    clear yNeg;
    clear yPos;
end

obj = obj1;
clear obj1;

for i = 1:r
    
    yNeg = zeros(c,d);
    yPos = zeros(c,d);
    
    for n = 3:c
        for d1 = 1:d
            yNeg(n, d1) = s_anti * exp(-alph) * obj(i,n-1,d1) + 2 * exp(-alph) * yNeg(n-1, d1) - exp(-2 * alph) * yNeg(n-2, d1) ;
        end
    end
    
    for n = c-2:-1:1
        for d1 = 1:d
            yPos(n, d1) = s_anti * exp(-alph) * obj(i,n+1,d1) + 2 * exp(-alph) * yPos(n+1, d1) - exp(-2 * alph) * yPos(n+2, d1) ;
        end
    end
    
    for n = 1:c
        for d1 = 1:d
            obj1(i,n,d1) = yNeg(n,d1) - yPos(n,d1);
        end
    end
    
end

Ix = obj1;
obj = obj1;
clear obj1;

gradImg = zeros(r,c);

for i=1:r
    for j=1:c
        for d1 = 1:d
            k = sqrt(Ix(i,j,d1)*Ix(i,j,d1) + Iy(i,j,d1)*Iy(i,j,d1));
            if(k > thresh)
                %checking for maxima along eta direction
                theta = atan(Iy(i,j,d1)/Ix(i,j,d1));
                eta = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                
                pt1 = [0;1];
                pt1 = eta*pt1;
                pt2 = [0;-1];
                pt2 = eta*pt2;
                
                lefti = round(i + pt1(1,1));
                leftj = round(j + pt1(2,1));
                righti = round(i + pt2(1,1));
                rightj = round(j + pt2(2,1));
                
                if(lefti > 0 && lefti <= r && righti > 0 && righti <= r && leftj > 0 && leftj <= c && rightj > 0 && rightj <= c)
                    if (sqrt(Ix(lefti,leftj,d1)*Ix(lefti,leftj,d1) + Iy(lefti,leftj,d1)*Iy(lefti,leftj,d1))<k && sqrt(Ix(righti,rightj,d1)*Ix(righti,rightj,d1) + Iy(righti,rightj,d1)*Iy(righti,rightj,d1))<k)
                        gradImg(i,j) = 255;
                    end
                end
            end
        end
    end
end

imshow(gradImg);
fprintf('Running time: %3.2d sec\n',cputime-timeTaken);
%imwrite(gradImg, strcat(strcat(file1(1:strfind(file1,'.')),'_edges3'),file1(strfind(file1,'.'):length(file1))));