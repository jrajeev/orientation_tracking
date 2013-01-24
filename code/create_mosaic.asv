function mosaic=create_mosaic(q_ukf)
    idx=find(tsimu_vicon);
    thresh_ang=50*pi/180;
    mosaic=zeros(2000,3000,3);
    tic
    for i=1:numel(idx)
        [angy,angz,angx]=quat2angle(transpose(q_fin(:,idx(i))));
        angy=angy*180/pi;
        angz=-angz;
        angy=-angy;
        %angx=-angx;
        if (abs(angy)<thresh_ang)    
            I=cam(:,:,:,tsimu_cam(idx(i)));
            I=double(imrotate(I,angx*180/pi));
            shiftr=10*tan(angy);
            shiftc=10*tan(angz);
            shiftr_std=1000-floor(size(I,1)/2);
            shiftc_std=1500-floor(size(I,2)/2);
            sr=floor(shiftr_std+shiftr);
            sc=floor(shiftc_std+shiftc);
            wtmatI=0.5*ones(size(I,1),size(I,2),3);
            wtmatI=(I>0).*wtmatI;
            mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)=uint8(((1-wtmatI).*mosaic(sr:sr+size(I,1)-1,sc:sc+size(I,2)-1,:)+wtmatI.*I(:,:,:)));
        end
    end
    G=fspecial('gaussian');
    mosaic(:,:,1)=conv2(mosaic(:,:,1),G,'same');
    mosaic(:,:,2)=conv2(mosaic(:,:,2),G,'same');
    mosaic(:,:,3)=conv2(mosaic(:,:,3),G,'same');
    img_idx=find(mosaic);
    [row,col,ch]=ind2sub(size(mosaic),img_idx);
    xxmin=min(row); xxmax=max(row);
    yymin=min(col); yymax=max(col);
    zzmin=min(ch); zzmax=max(ch);
    toc
    figure,imshow(uint8(mosaic(xxmin:xxmax,yymin:yymax,:)));
end