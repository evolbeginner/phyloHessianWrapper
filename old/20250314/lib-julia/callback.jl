function callback(x)
	global prev_value, prev_g_norm

    for field in fieldnames(typeof(x))
        #println(field,"\t",getfield(x, field))
		;
    end 

    curr_value = x.value
	curr_g_norm = x.g_norm

    convergence_criterion = abs(curr_g_norm - prev_g_norm) < 1e-2
	#println("diff","\t",abs(curr_g_norm - prev_g_norm))
    #convergence_criterion = abs(curr_value - prev_value) < 1e-2
	println("diff","\t",abs(curr_value - prev_value))

	#inv_hessian_str = match(r"~inv\(H\): (.+)", string(x)).captures[1]
	#println( inv(eval(Meta.parse(matrix_string))) )

    println()

    if convergence_criterion
        println("Optimization converged!")
        return true  # Returning true stops the optimization
    end 

    # Update previous result for the next iteration
	prev_value = curr_value
    prev_g_norm = curr_g_norm
    return false
end


