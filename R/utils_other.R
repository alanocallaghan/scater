.switch_arg_names <- function(old.val, new.val) {
   if (!is.null(old.val)) {
        old.arg <- deparse(substitute(old.val))
        new.arg <- deparse(substitute(new.val))
        .Deprecated(new=new.arg, old=old.arg)
        old.val
    } else {
        new.val
   }
}
