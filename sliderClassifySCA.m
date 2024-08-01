function h = sliderClassifySCA(RGB)

%determine RGB colorspace (only 8-bit or 16-bit rgb right now)
if isa(RGB,"uint8")
    M = 2^8 - 1;
elseif isa(RGB,"uint16")
    M = 2^16 - 1;
else
    error('Unsupported image datatype')
end

%prepare guidata struct 
h=struct;
h.f=figure("Position",[100 100 900 900]);clf(h.f)%open a clean figure
h.ax=axes('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.2 0.8 0.75]);

%prepare data
h.RGBc = RGB;
h.SCA = zeros(size(h.RGBc,[1 2]));

%compute metrics
R = single(h.RGBc(:,:,1))./M; G = single(h.RGBc(:,:,2))./M; B = single(h.RGBc(:,:,3))./M;
%whiteness
W = (((1 - R).^2) + ((1 - G).^2) + ((1 - B).^2)).^0.5;
W = W./((3*(1.^2)).^0.5); %normalize based on max possible distance
W = 1 - W; %invert so high values indicate whiter areas
h.W = W;
%blue - red band difference
h.BRd = B - R;
%get minimum and maximum values for setting slider bounds
v1 = min([h.W(:); h.BRd(:)]);
v2 = max([h.W(:); h.BRd(:)]);

%initialize plot and mask
h.im = imagesc(h.RGBc,'Parent',h.ax);
h.alphamask = alphamask(h.SCA, [0 0 0],0.75,h.ax);
title(h.ax,{"(1) Adjust sliders until only snowy areas are masked",...
    "(2) Press continue to store as a reference SCA map"},'FontSize',16);

%create sliders
%slider to control weighting parameter
h.slider_w = uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.05 0.8 0.09],...
    'Style','slider',...
    'Min',0,'Max',1,'Value',0.5,...
    'Callback',@slider_callback);
h.slider_label1 = annotation('textbox',[0.05 0.1 0.05 0.05],...
    'String','w (low)','Edgecolor','none','FontWeight','bold');
h.slider_label2 = annotation('textbox',[0.90 0.1 0.05 0.05],...
    'String','w (high)','Edgecolor','none','FontWeight','bold');

%slider to control blue threshold
h.slider_SCM = uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.1 0.01 0.8 0.09],...
    'Style','slider',...
    'Min',v1,'Max',v2,'Value',(v1+v2)/2,...
    'Callback',@slider_callback);
h.slider_label3 = annotation('textbox',[0.03 0.05 0.05 0.05],...
    'String','SCMthresh (low)','Edgecolor','none','FontWeight','bold');
h.slider_label4 = annotation('textbox',[0.90 0.05 0.05 0.05],...
    'String','SCMthresh (high)','Edgecolor','none','FontWeight','bold');

%create button to save settings when finalized
h.btn = uicontrol('Parent',h.f,...
    'Units','Normalized',...
    'Position',[0.45 0.03 0.10 0.05],...
    'Style','pushbutton',...
    'FontWeight','bold',...
    'String','Continue',...
    'Callback',@button_fcn);

%store handles and data to guidata
guidata(h.f,h)
waitfor(h.f)

    function slider_callback(hObject,eventdata)

        %retrieve struct
        h=guidata(hObject);

        %get all slider values
        w = h.slider_w.Value;
        SCM_thresh = h.slider_SCM.Value;

        %compute snow cover metric using the weighted combination of W and BRd
        SCM = h.W .*(w) + h.BRd .* (1 - w);

        %update SCA map using combined metrics
        %set(h.SCA,single(SCM > SCM_thresh));
        h.SCA = single(SCM > SCM_thresh);

        %update figure with new SCA mask
        set(h.alphamask,"AlphaData",h.SCA * 0.75);
        
        %update guidata
        guidata(h.f,h)

    end

    function button_fcn(hObject,eventdata)

        %close figure and return cleaned structure
        delete(gcf);
        h = rmfield(h,'f');
        h = rmfield(h,'ax');
        h = rmfield(h,'im');
        h = rmfield(h,'alphamask');
        h = rmfield(h,'slider_w');
        h = rmfield(h,'slider_label1');
        h = rmfield(h,'slider_label2');
        h = rmfield(h,'slider_label3');
        h = rmfield(h,'slider_label4');
        h = rmfield(h,'slider_SCM');
        h = rmfield(h,'btn');
        return

    end

end
