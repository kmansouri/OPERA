function pred = pack_pred_struct(D, neighbors, w, dc, varargin)

pred = struct();
pred.D = D;
pred.neighbors = neighbors;
pred.w = w;
pred.dc = dc;

if nargin >= 5 && ~isempty(varargin{1})
    pred.class_pred = varargin{1};
end

if nargin >= 6 && ~isempty(varargin{2})
    pred.y_pred_weighted = varargin{2};
end

end