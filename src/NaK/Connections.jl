module Connections

export chain, kchain, ring
export sheet, torus

function chain(n)
	connections = zeros(n,n)
	for i = 1:n
		if i == 1
			connections[1,2] = 1
		elseif i == n
			connections[n-1,n] = 1
		else
			connections[i,i-1] = 1
			connections[i,i+1] = 1
		end
	end
	return connections
end

function kchain(n,k)
	connections = zeros(n,n)
	for i = 1:n
		for j = (i-k):(i+k)
			if i != j
				if j >= 1 && j <= n
					connections[i,j] = 1/(abs(i-j)^2)
				end
			end
		end
	end
	return connections
end

function ring(n)
	connections = zeros(n,n)
	for i = 1:n
		if i == 1
			j = n
			k = 2
		elseif i == n
			j = n - 1
			k = 1
		else
			j = i - 1
			k = i + 1
		end
		connections[i,j] = 1
		connections[i,k] = 1
	end
	return connections
end


function sheet(n)
	sqrtn = Int64(floor(sqrt(n)))
	connections = zeros(n,n)
	for i = 1:n
		if i == 1
			connections[i,i+1] = 1
		elseif i == n
			connections[i,i-1] = 1
		else
			connections[i,i+1] = 1
			connections[i,i-1] = 1
		end
		u = i - sqrtn
		if u > 0
			connections[i,u] = 1
		end
		d = i + sqrtn
		if d <= n
			connections[i,d] = 1
		end
	end
	return connections
end

function torus(n)
	sqrtn = Int64(floor(sqrt(n)))
	connections = zeros(n,n)
	for i = 1:n
		if i == 1
			l = n
			r = 2
		elseif i == n
			l = n - 1
			r = 1
		else
			l = i - 1
			r = i + 1
		end
		u = i - sqrtn
		if u <= 0
			u = n + u
		end
		d = i + sqrtn
		if d > n
			d = d - n
		end

		connections[i,u] = 1
		connections[i,d] = 1
		connections[i,l] = 1
		connections[i,r] = 1
	end
	return connections
end

end
