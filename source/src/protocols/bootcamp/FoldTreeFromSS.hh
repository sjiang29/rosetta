class FoldTreeFromSS {
public:
    FoldTreeFromSS( std::string const & ssstring );

    core::kinematics::FoldTree const &
    fold_tree() const;

    protocols::loops::Loop const &
    loop( core::Size index ) const;

    core::Size
    loop_for_residue( core::Size seqpos ) const;

private:
    core::kinematics::FoldTree ft_;
    utility::vector1< protocols::loops::Loop > loop_vector_;
    utility::vector1< core::Size > loop_for_residue_;
};
