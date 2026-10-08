# Shared task helpers: run a threaded loop so that an error in one task
# reaches the caller as the exception the task itself raised.
#
# `Threads.@threads`, `@sync` and `wait` wrap a task's error in a
# `TaskFailedException` (inside a `CompositeException` for `@sync`), so a
# caller writing `catch e; e isa ArgumentError` sees it on a serial path and
# misses it on a threaded one. The parallel samplers of ERGM.jl and
# ERGMRank.jl and the bootstrap loop here all need the same two lines; they
# live here once.

"""
    unwrap_task_exception(e) -> Exception

The first underlying exception of a `TaskFailedException` /
`CompositeException` nest (what `wait`, `fetch`, `@sync` and
`Threads.@threads` raise when a task fails), or `e` itself when it is
neither. Peels any depth of nesting, taking the first exception of a
`CompositeException`.

# Example
```julia
import NetworkCore: unwrap_task_exception
t = Threads.@spawn throw(ArgumentError("bad input"))
e = try wait(t) catch err; err end
e isa TaskFailedException                        # true
unwrap_task_exception(e)                         # ArgumentError("bad input")
unwrap_task_exception(ErrorException("x"))       # returned unchanged
```
"""
function unwrap_task_exception(e)
    while true
        if e isa TaskFailedException
            e = e.task.exception
        elseif e isa CompositeException && !isempty(e.exceptions)
            e = first(e.exceptions)
        else
            return e
        end
    end
end

"""
    spawn_all(f, n::Integer) -> nothing

Run `f(1), f(2), …, f(n)` on separate tasks (`Threads.@spawn`) and wait for
**all** of them. If any task fails, the first failure in task order is
rethrown as the task's own exception (through
[`unwrap_task_exception`](@ref)), never as a `TaskFailedException`, and only
after every task has finished, so none is left running behind the error.

The failing task's backtrace is not lost: the `TaskFailedException` is kept on
the exception stack as the cause of the rethrown error. In a `catch` block,
`current_exceptions()` returns both, and the first entry's
`task.backtrace` holds the frames inside the task; an uncaught error prints
them under "caused by: TaskFailedException … nested task error".

The order in which tasks run is up to the scheduler, so `f(c)` must write
only to its own slot of any shared output (`out[c] = …`). Seeding is the
caller's job: draw one seed per task from the caller's `rng` before calling,
so the result does not depend on the thread count.

# Example
```julia
import NetworkCore: spawn_all
out = zeros(4)
spawn_all(c -> (out[c] = c^2), 4)
out                                                         # [1.0, 4.0, 9.0, 16.0]
try
    spawn_all(c -> c == 3 ? throw(ArgumentError("task 3")) : nothing, 4)
catch e
    e                                                       # ArgumentError("task 3")
    cause = first(current_exceptions()).exception           # the TaskFailedException
    cause.task.backtrace                                    # frames inside task 3
end
```
"""
function spawn_all(f::F, n::Integer) where {F}
    tasks = [Threads.@spawn f(c) for c in 1:n]
    err = nothing
    for t in tasks
        try
            wait(t)
        catch e
            err === nothing && (err = e)
        end
    end
    err === nothing && return nothing
    # Raise the task's own exception from inside a handler for the wrapper, so
    # the wrapper stays on the exception stack as its cause: the
    # `TaskFailedException` carries the failing task and the backtrace from
    # inside it, which `current_exceptions()` returns and an uncaught error
    # prints under "caused by: … nested task error". Thrown from outside a
    # handler, the error's only backtrace would start here, in `spawn_all`.
    try
        throw(err)
    catch
        inner = unwrap_task_exception(err)
        inner === err ? rethrow() : throw(inner)
    end
end
