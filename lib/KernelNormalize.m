function K2 = KernelNormalize(K1, type)
if ~exist('type', 'var') || isempty(type)
    type = '';
end

switch lower(type)
    case lower('')
        K2 = K1;
    case lower('NCW')
        nSmp = size(K1, 1);
        D_mhalf = sum(K1,2).^-.5;
        D_mhalf = spdiags(D_mhalf,0,nSmp,nSmp);
        K2 = D_mhalf * K1 * D_mhalf;
		clear K1;
	case lower('MAX')
		K2 = K1 / max(K1(:));
		clear K1;
	case lower('uni_trace')
		K2 = K1 / trace(K1);
		clear K1;
	case lower('Sym')
		K2 = (K1 + K1')/2;
		K2 = real(K2);
		clear K1;
	case lower('SCALE')
		min_k = min(K1(:));
		max_k = max(K1(:));
		diff = max(max_k - min_k, eps);
		K2 = (K1 - min_k) / diff;
		clear K1;
	case lower('NCW-SCALE-SYM')
		K11 = KernelNormalize(K1, 'NCW');
		K12 = KernelNormalize(K11, 'SCALE');
		K2 = KernelNormalize(K12, 'SYM');
		clear K11 K12 K1;
    case lower('Sample-Scale');
        kw = max(diag(K1), eps);
        kw = kw .^-0.5;
        K2 = bsxfun(@times, K1, kw);
        K2 = bsxfun(@times, K2, kw');
end
end