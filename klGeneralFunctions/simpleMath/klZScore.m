function out = klZScore(in,params)

if size(params,1) == 1,
    params = repmat(params,size(in,1),1);
elseif size(params,1) ~= size(in,1),
    error('Inputs must have equal lengths');
end

params(params(:,2) == 0) = nan;

out = (in-repmat(params(:,1),1,size(in,2))).*repmat(1./params(:,2),1,size(in,2));
% for im = 1:size(in,1),
%     if mod(im,1000) == 0,
%         fprintf('%d of %d\n',im,size(in,1));
%     end
%     out(im,:) = (in(im,:)-params(im,1))./params(im,2);
%     if any(out(im,:) == -inf),
%         keyboard
%     end
% end
