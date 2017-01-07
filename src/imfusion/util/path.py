def symlink_relative(src_path, dest_path):
    dest_path.symlink_to(src_path.relative_to(dest_path.parent))
